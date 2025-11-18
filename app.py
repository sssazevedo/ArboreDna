import os
import re
import unicodedata
import string
import pandas as pd
import networkx as nx
from collections import deque
from flask import Flask, render_template, request
import json
from networkx.algorithms.community import greedy_modularity_communities

# ------------------------------------------------------------
# CONFIGURAÇÕES DO FLASK
# ------------------------------------------------------------
app = Flask(__name__)
app.secret_key = 'arbore_dna_v1'
UPLOAD_FOLDER = "uploads"
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

# ------------------------------------------------------------
# FUNÇÕES DE NORMALIZAÇÃO
# ------------------------------------------------------------
def norm_name(s: str) -> str:
    if s is None:
        return ""
    s = str(s).strip().replace("\n", " ").replace("\r", " ")
    s = "".join(ch for ch in s if ch.isprintable())
    if s == "":
        return ""
    original = s
    s = unicodedata.normalize("NFKC", s)
    s = re.sub(r"[^0-9A-Za-zÁ-ú' -]", " ", s)
    s = " ".join(s.split())
    parts = s.split()
    if len(parts) == 1:
        return parts[0].lower()
    s = unicodedata.normalize("NFKD", s)
    s = "".join(c for c in s if not unicodedata.combining(c))
    s = "".join(ch if (ch.isalnum() or ch in " -'") else " " for ch in s)
    s = " ".join(s.split()).lower()
    return s

# ------------------------------------------------------------
# PARSE DE SEGMENTOS (PARA USUÁRIO E PAIS)
# ------------------------------------------------------------
def parse_segments_csv(path):
    try:
        df = pd.read_csv(path, encoding='utf-8', skipinitialspace=True)
    except:
        df = pd.read_csv(path, encoding='latin-1', skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    name_col = next((c for c in df.columns if c.lower() in ("name","matchedname","match","kit","match name")), None)
    chr_col  = next((c for c in df.columns if c.lower() in ("chr","chrom","chromosome")), None)
    s_col    = next((c for c in df.columns if c.lower() in ("start","startpos","start position","b37start")), None)
    e_col    = next((c for c in df.columns if c.lower() in ("end","endpos","end position","b37end")), None)
    cm_col   = next((c for c in df.columns if "cm" in c.lower()), None)
    if not (name_col and chr_col and s_col and e_col and cm_col):
        raise ValueError("Arquivo de segmentos inválido.")
    segs_by = {}
    for _, r in df.iterrows():
        nm = r.get(name_col, "")
        if pd.isna(nm): continue
        nm = str(nm).strip()
        if nm == "": continue
        try:
            chr_val = str(r[chr_col]).strip().upper()
            chrn = int(chr_val.replace("X", "23"))
        except: continue
        try:
            s = int(float(r[s_col]))
            e = int(float(r[e_col]))
            cm = float(r[cm_col])
        except: continue
        if cm < 7: continue
        segs_by.setdefault(nm, []).append({"chr": chrn, "start": s, "end": e, "cm": cm})
    return segs_by

def parse_parent_segment_files(file_list):
    """Parses segment files and returns a simple set of normalized match names."""
    names_set = set()
    for f in file_list:
        if f.filename == "": continue
        p = os.path.join(UPLOAD_FOLDER, f.filename)
        f.save(p)
        try:
            parsed = parse_segments_csv(p)
            for name in parsed.keys():
                nn = norm_name(name)
                if nn:
                    names_set.add(nn)
        except Exception as e:
            print(f"[AVISO] Falha ao ler arquivo de segmento do pai/mãe: {e}")
    return names_set

# ------------------------------------------------------------
# PARSE DE TRIANGULAÇÃO (COM PESOS)
# ------------------------------------------------------------
def parse_triang_csv(path):
    try:
        df = pd.read_csv(path, encoding='utf-8', skipinitialspace=True)
    except:
        df = pd.read_csv(path, encoding='latin-1', skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    a_col = next((c for c in df.columns if any(k in c.lower() for k in ("match a","match1","name a","person a","kit1 name"))), None)
    b_col = next((c for c in df.columns if any(k in c.lower() for k in ("match b","match2","name b","person b","kit2 name"))), None)
    chr_col = next((c for c in df.columns if c.lower() in ("chr","chrom","chromosome")), None)
    cm_col   = next((c for c in df.columns if "cm" in c.lower()), None)
    if not (a_col and b_col and chr_col and cm_col):
        print("[AVISO] Colunas de Triangulação (A, B, Chr, cM) não encontradas. Pulando arquivo.")
        return {}
    tri_with_cm = {}
    for _, r in df.iterrows():
        a = norm_name(str(r[a_col]))
        b = norm_name(str(r[b_col]))
        try:
            chr_val = str(r[chr_col]).strip().upper()
            chrn = int(chr_val.replace("X", "23"))
        except: continue
        try:
            cm = float(r[cm_col])
        except: cm = 0.0
        if not a or not b or a == b or cm < 7.0: continue
        if a > b: a, b = b, a
        tkey = (a, b, chrn)
        tri_with_cm[tkey] = max(tri_with_cm.get(tkey, 0.0), cm)
    print(f"[DEBUG] 'parse_triang_csv' finalizou. Pares ponderados lidos: {len(tri_with_cm)}")
    return tri_with_cm

# ------------------------------------------------------------
# OVERLAP
# ------------------------------------------------------------
def overlap(seg1, seg2, min_cm=7.0):
    if seg1["chr"] != seg2["chr"]: return False
    if seg1["end"] < seg2["start"] or seg2["end"] < seg1["start"]: return False
    return seg1["cm"] >= min_cm and seg2["cm"] >= min_cm

# ------------------------------------------------------------
# CONSTRUIR GRAFO (COM PESOS)
# ------------------------------------------------------------
def build_side_graph(match_names, segs_by, tri_pairs_with_cm=None):
    tri_pairs_with_cm = tri_pairs_with_cm or {}
    tri_pairs_provided = len(tri_pairs_with_cm) > 0
    
    mset = {norm_name(n) for n in match_names}
    by_chr = {}
    for name, segs in segs_by.items():
        nn = norm_name(name)
        if nn not in mset: continue
        for s in segs:
            by_chr.setdefault(s["chr"], {}).setdefault(nn, []).append(s)
            
    out = {}
    
    for chrn, mp in by_chr.items():
        same_weighted = {}
        opp  = set()
        names = sorted(mp.keys())
        
        for i in range(len(names)):
            for j in range(i+1, len(names)):
                a, b = names[i], names[j]
                
                has_overlap = False
                for sa in mp[a]:
                    for sb in mp[b]:
                        if overlap(sa, sb):
                            has_overlap = True
                if not has_overlap: continue
                
                tkey_chr = (a,b,chrn) if a<b else (b,a,chrn)
                tkey = (a,b) if a<b else (b,a)

                # --- LÓGICA CORRIGIDA ---
                if tri_pairs_provided:
                    # Se temos triangulação, SÓ usamos links triangulados.
                    if tkey_chr in tri_pairs_with_cm:
                        cm_weight = tri_pairs_with_cm[tkey_chr]
                        # Acumula o peso (soma do QUADRADO dos cMs)
                        same_weighted[tkey] = same_weighted.get(tkey, 0.0) + (cm_weight * cm_weight)
                    
                    # O 'else' (de peso 1.0) FOI REMOVIDO para não poluir o grafo.
                
                else:
                    # Se NÃO temos triangulação, usamos o overlap como link (peso 1.0)
                    same_weighted[tkey] = same_weighted.get(tkey, 0.0) + 1.0
        
        if same_weighted or opp:
            out[chrn] = {"same": same_weighted, "opp": opp}
            
    return out

# ------------------------------------------------------------
# FUNÇÃO AUXILIAR PARA CONSTRUIR ÁRVORE (MODO HIERÁRQUICO - CORRIGIDO)
# ------------------------------------------------------------
def build_tree_json(phased_groups, kept, G, use_weight=True, user_resolution=1.0):
    tree = {"nodes": [], "edges": []}
    table_groups = [] # <-- NOVO: Lista para a tabela
    
    # Define o parâmetro de peso para o Louvain
    louvain_weight = 'weight' if use_weight else None
    
    colors = ["#F44336", "#2196F3", "#4CAF50", "#FFC107", "#9C27B0", "#00BCD4"]
    tree["nodes"].append({"id": "root", "label": "Você", "side": None})
    
    main_branches = [
        ("paternal", "Lado Paterno", 1), # Azul
        ("maternal", "Lado Materno", 0), # Vermelho
        ("both", "Ambos os Lados", 2), # Verde
        ("unknown", "Desconhecido", 3)  # Amarelo
    ]
    
    # --- MODO FALLBACK ---
    if phased_groups is None:
        # ------------------------------------------------------------
        # CORREÇÃO: Criar uma cópia e remover isolados ANTES de clusterizar
        # ------------------------------------------------------------
        G_sub_limpo = G.copy() # Criar cópia para não modificar o 'G' original
        nodes_isolados = list(nx.isolates(G_sub_limpo))
        if nodes_isolados:
            G_sub_limpo.remove_nodes_from(nodes_isolados)
            print(f"[INFO] build_tree (Fallback): Removendo {len(nodes_isolados)} isolados.")
        # ------------------------------------------------------------
        
        G_sub = G_sub_limpo # Usar o grafo limpo
        subramos = list(greedy_modularity_communities(G_sub, resolution=0.5, weight=louvain_weight))
        
        if not subramos and list(G_sub.nodes()): # Fallback
            subramos = [set(G_sub.nodes())]

        for i, subramo_nodes in enumerate(subramos):
            subramo_id = f"subramo_{i}"
            subramo_cor = i % len(colors)
            subramo_label = f"Subramo {i}"
            
            tree["nodes"].append({"id": subramo_id, "label": subramo_label, "side": subramo_cor})
            tree["edges"].append({"from": "root", "to": subramo_id})

            # 2. Encontrar Ramos
            G_ramo = G_sub.subgraph(subramo_nodes)
            ramos = list(greedy_modularity_communities(G_ramo, resolution=user_resolution, weight=louvain_weight))

            if not ramos and list(G_ramo.nodes()): # Fallback
                 ramos = [set(G_ramo.nodes())]
            
            if len(ramos) == 1:
                ramo_nodes = list(ramos)[0]
                ramo_id_para_link = subramo_id
                ramo_label = subramo_label # O nome do grupo é o do Subramo
                ramo_nodes_flat = True
            else:
                ramo_nodes = None
                ramo_nodes_flat = False

            for j, ramo_nodes_loop in enumerate(ramos):
                current_matches = [] # <-- NOVO: Lista de matches para este grupo
                
                if not ramo_nodes_flat:
                    ramo_nodes_inner = ramo_nodes_loop
                    ramo_id = f"ramo_{i}_{j}"
                    ramo_label = f"Ramo {i}-{j}" # Nome do grupo é o Ramo
                    tree["nodes"].append({"id": ramo_id, "label": ramo_label, "side": subramo_cor})
                    tree["edges"].append({"from": subramo_id, "to": ramo_id})
                    ramo_id_para_link = ramo_id
                else:
                    ramo_nodes_inner = ramo_nodes
                
                for member in ramo_nodes_inner:
                    if member not in kept: continue
                    orig_name, cmval = kept[member]
                    tree["nodes"].append({
                        "id": member,
                        "label": f"{orig_name} ({cmval:.1f} cM)",
                        "side": subramo_cor,
                        "cm": cmval
                    })
                    tree["edges"].append({"from": ramo_id_para_link, "to": member})
                    current_matches.append({"name": orig_name, "cm": cmval}) # <-- NOVO
                
                # Adiciona o grupo e seus matches à tabela
                if current_matches:
                    table_groups.append({"group_name": ramo_label, "matches": current_matches})

                if ramo_nodes_flat:
                    break
        
        return tree, table_groups # <-- MUDANÇA: Retorna a tabela

    # --- MODO PHASED ---
    for group_key, group_label, group_color_index in main_branches:
        
        group_nodes = phased_groups[group_key]
        if not group_nodes: continue
        
        branch_id = f"branch_{group_key}"
        link_from = "root"
        
        tree["nodes"].append({"id": branch_id, "label": group_label, "side": group_color_index})
        tree["edges"].append({"from": link_from, "to": branch_id})
        link_from = branch_id

        # ------------------------------------------------------------
        # CORREÇÃO: Criar um subgrafo-cópia e remover isolados ANTES
        # ------------------------------------------------------------
        # Criamos uma CÓPIA do subgrafo, não uma view, para poder modificá-lo
        G_sub_limpo = G.subgraph(group_nodes).copy() 
        
        nodes_isolados_no_subgrafo = list(nx.isolates(G_sub_limpo))
        if nodes_isolados_no_subgrafo:
            G_sub_limpo.remove_nodes_from(nodes_isolados_no_subgrafo)
            print(f"[INFO] build_tree ({group_key}): Removendo {len(nodes_isolados_no_subgrafo)} isolados.")
        # ------------------------------------------------------------

        G_sub = G_sub_limpo # Usar o grafo limpo
        subramos = list(greedy_modularity_communities(G_sub, resolution=0.5, weight=louvain_weight))
        
        if not subramos and list(G_sub.nodes()):
            subramos = [set(G_sub.nodes())]

        for i, subramo_nodes in enumerate(subramos):
            
            subramo_id = f"subramo_{group_key}_{i}"
            subramo_label = f"{group_label} - Subramo {i}"
            tree["nodes"].append({"id": subramo_id, "label": f"Subramo {i}", "side": group_color_index})
            tree["edges"].append({"from": link_from, "to": subramo_id})

            G_ramo = G_sub.subgraph(subramo_nodes)
            ramos = list(greedy_modularity_communities(G_ramo, resolution=user_resolution, weight=louvain_weight))

            if not ramos and list(G_ramo.nodes()):
                 ramos = [set(G_ramo.nodes())]
            
            if len(ramos) == 1:
                ramo_nodes = list(ramos)[0]
                ramo_id_para_link = subramo_id
                ramo_label = subramo_label # Nome do grupo é o Subramo
                ramo_nodes_flat = True
            else:
                ramo_nodes = None
                ramo_nodes_flat = False

            for j, ramo_nodes_loop in enumerate(ramos):
                current_matches = [] # <-- NOVO
                
                if not ramo_nodes_flat:
                    ramo_nodes_inner = ramo_nodes_loop
                    ramo_id = f"ramo_{group_key}_{i}_{j}"
                    ramo_label = f"{subramo_label} - Ramo {j}" # Nome do grupo
                    tree["nodes"].append({"id": ramo_id, "label": f"Ramo {i}-{j}", "side": group_color_index})
                    tree["edges"].append({"from": subramo_id, "to": ramo_id})
                    ramo_id_para_link = ramo_id
                else:
                    ramo_nodes_inner = ramo_nodes
                
                for member in ramo_nodes_inner:
                    if member not in kept: continue
                    orig_name, cmval = kept[member]
                    tree["nodes"].append({
                        "id": member,
                        "label": f"{orig_name} ({cmval:.1f} cM)",
                        "side": group_color_index,
                        "cm": cmval
                    })
                    tree["edges"].append({"from": ramo_id_para_link, "to": member})
                    current_matches.append({"name": orig_name, "cm": cmval}) # <-- NOVO
                
                if current_matches:
                    table_groups.append({"group_name": ramo_label, "matches": current_matches})

                if ramo_nodes_flat:
                    break
        
    return tree, table_groups
# ------------------------------------------------------------
# ROTA PRINCIPAL
# ------------------------------------------------------------
@app.route("/", methods=["GET","POST"])
def index():
    if request.method == "POST":

        # ===============================
        # 1) CARREGAR ARQUIVOS
        # ===============================
        seg_files_user = request.files.getlist("segment_files")
        seg_files_pai = request.files.getlist("father_segment_files")
        seg_files_mae = request.files.getlist("mother_segment_files")
        tri_files = request.files.getlist("triangulation_csv")

        if not seg_files_user or seg_files_user[0].filename == "":
            return render_template("index.html", error="Envie ao menos 1 CSV de segmentos do Usuário.")
        
        is_phased = (seg_files_pai and seg_files_pai[0].filename != "") and \
                    (seg_files_mae and seg_files_mae[0].filename != "")
                    
        # --- NOVO: Capturar Resolução do Formulário ---
        try:
            # Pega o valor do <select>, usa '1.0' como padrão se falhar
            user_resolution = float(request.form.get("resolution_select", "1.0"))
        except ValueError:
            user_resolution = 1.0 # Padrão de segurança
        print(f"[INFO] Usando Resolução de Ramos: {user_resolution}")
        # ------------------------------------------------

        # ===============================
        # 2) PARSEAR SEGMENTOS DO USUÁRIO
        # ===============================
        raw_segs = {}
        for f in seg_files_user:
            if f.filename == "": continue
            p = os.path.join(UPLOAD_FOLDER, f.filename)
            f.save(p)
            parsed = parse_segments_csv(p)
            for name, segs in parsed.items():
                raw_segs.setdefault(name, []).extend(segs)

        # ===============================
        # 3) PARSEAR TRIANGULAÇÃO (SEMPRE NECESSÁRIO)
        # ===============================
        all_tri = {}
        if tri_files and tri_files[0].filename != "":
            for f in tri_files:
                if f.filename == "": continue
                p = os.path.join(UPLOAD_FOLDER, f.filename)
                f.save(p)
                parsed_tri = parse_triang_csv(p)
                all_tri.update(parsed_tri)
                
        tem_triangulacao = len(all_tri) > 0
        
        # ===============================
        # 4) NORMALIZAR E FILTRAR USUÁRIO
        # ===============================
        normalized_segs = {}
        norm_to_orig = {}
        for orig_name, segs in raw_segs.items():
            nn = norm_name(orig_name)
            if not nn: continue
            normalized_segs.setdefault(nn, []).extend(segs)
            if nn not in norm_to_orig:
                norm_to_orig[nn] = orig_name
        
        if not normalized_segs:
            return render_template("index.html", error="Nenhum segmento válido após normalização.")

        total_cm = {}
        for norm_name_key, segs in normalized_segs.items():
            total = sum(s["cm"] for s in segs)
            total_cm[norm_name_key] = (norm_to_orig.get(norm_name_key, norm_name_key), total)

        CM_THRESHOLD = 30
        kept = {n:(orig,cm) for n,(orig,cm) in total_cm.items() if cm >= CM_THRESHOLD}
        if not kept:
            return render_template("index.html", error=f"Nenhum match >= {CM_THRESHOLD} cM.")
        
        normalized_segs = {n: normalized_segs[n] for n in kept}
        names = list(normalized_segs.keys())

        # ===============================
        # 5) CONSTRUIR GRAFO PONDERADO (SEMPRE NECESSÁRIO)
        # ===============================
        sg = build_side_graph(names, normalized_segs, all_tri)
        same_all_weighted = {}
        nodes = set()
        for chrn, edges in sg.items():
            for (a, b), weight in edges["same"].items():
                if a in kept and b in kept:
                    same_all_weighted[(a,b)] = same_all_weighted.get((a,b), 0.0) + weight
                    nodes.add(a); nodes.add(b)

        G = nx.Graph()
        G.add_nodes_from(kept)
        
        # Vamos verificar se temos dados de triangulação reais
        if tem_triangulacao:
            print("[INFO] Grafo Ponderado (Triangulação) ATIVADO.")
            for (a, b), weight in same_all_weighted.items():
                if a in kept and b in kept:
                    G.add_edge(a, b, weight=weight)
        else:
            # Modo Fallback (Sem Triangulação): Grafo não-ponderado
            print("[INFO] Grafo NÃO Ponderado (Overlap) ATIVADO.")
            for (a, b) in same_all_weighted.keys(): # Ignora o 'weight=1.0'
                if a in kept and b in kept:
                    G.add_edge(a, b) # Adiciona a aresta SEM peso
        
               
        # ===============================
        # 6) ESCOLHER O MODO E GERAR ÁRVORE
        # ===============================
        phased_groups = None
        tree = {}
        
        if is_phased:
            print("[INFO] Modo Phased Ativado.")
            matches_do_pai = parse_parent_segment_files(seg_files_pai)
            matches_da_mae = parse_parent_segment_files(seg_files_mae)

            # --- ETAPA 6A: Seeds ---
            initial_groups = {
                "paternal": set(),
                "maternal": set(),
                "both": set(),
                "unknown": set()
            }

            for nn_key in kept.keys():
                in_pai = nn_key in matches_do_pai
                in_mae = nn_key in matches_da_mae
                if in_pai and in_mae:
                    initial_groups["both"].add(nn_key)
                elif in_pai:
                    initial_groups["paternal"].add(nn_key)
                elif in_mae:
                    initial_groups["maternal"].add(nn_key)
                else:
                    initial_groups["unknown"].add(nn_key)

            print(f"[INFO] Phasing (Inicial): {len(initial_groups['paternal'])} Pat / "
                  f"{len(initial_groups['maternal'])} Mat / "
                  f"{len(initial_groups['both'])} Ambos / "
                  f"{len(initial_groups['unknown'])} Desconhecidos")

            # --- ETAPA 6B: Phasing Final ---
            phased_groups = {
                "paternal": set(initial_groups["paternal"]),
                "maternal": set(initial_groups["maternal"]),
                "both": set(),
                "unknown": set()
            }

            matches_to_reclassify = initial_groups["unknown"].union(initial_groups["both"])

            for u_match in matches_to_reclassify:
                weight_to_pat = 0.0
                weight_to_mat = 0.0

                if u_match not in G:
                    phased_groups["unknown"].add(u_match)
                    continue

                for v_neighbor, edge_data in G[u_match].items():
                    w = edge_data.get("weight", 0.0)
                    if v_neighbor in initial_groups["paternal"]:
                        weight_to_pat += w
                    elif v_neighbor in initial_groups["maternal"]:
                        weight_to_mat += w

                if weight_to_pat > weight_to_mat:
                    phased_groups["paternal"].add(u_match)
                elif weight_to_mat > weight_to_pat:
                    phased_groups["maternal"].add(u_match)
                else:
                    if u_match in initial_groups["both"]:
                        phased_groups["both"].add(u_match)
                    else:
                        phased_groups["unknown"].add(u_match)

            print(f"[INFO] Phasing (Final): {len(phased_groups['paternal'])} Pat / "
                  f"{len(phased_groups['maternal'])} Mat / "
                  f"{len(phased_groups['both'])} Ambos / "
                  f"{len(phased_groups['unknown'])} Desconhecidos")

            # ------------------------------------------------------------
            # B2 — REFORÇO ENTRE MEMBROS DO MESMO LADO (AGORA NO LOCAL CERTO)
            # ------------------------------------------------------------
            print("[B2] Reforçando conexões internas entre membros do mesmo lado...")

            for lado in ["paternal", "maternal", "both"]:
                nodes_lado = list(phased_groups[lado])
                for i in range(len(nodes_lado)):
                    for j in range(i + 1, len(nodes_lado)):
                        a = nodes_lado[i]
                        b = nodes_lado[j]
                        if a in G and b in G and not G.has_edge(a, b):
                            G.add_edge(a, b, weight=1.0)

            # ------------------------------------------------------------
            # B3 — MATCHES FORTES (>50 cM)
            # ------------------------------------------------------------
            print("[B3] Aplicando reforço para matches fortes (>50 cM)...")

            cm_por_match = {n: kept[n][1] for n in kept}

            strong = [m for m in cm_por_match if cm_por_match[m] >= 50]

            for i in range(len(strong)):
                for j in range(i + 1, len(strong)):
                    a = strong[i]
                    b = strong[j]
                    if not G.has_edge(a, b):
                        G.add_edge(a, b, weight=0.5)

            # ------------------------------------------------------------
            # Construir árvore já com B2+B3 aplicados
            # ------------------------------------------------------------
            tree, table_data = build_tree_json(phased_groups, kept, G, use_weight=tem_triangulacao, user_resolution=user_resolution)


        else:
            # --- MODO FALLBACK (LOUVAIN PONDERADO) ---
            print("[INFO] Modo Fallback (Louvain Ponderado) Ativado.")
            
            # Passamos 'None' para os grupos, e a função irá rodar o Louvain em tudo
            tree, table_data = build_tree_json(None, kept, G, use_weight=tem_triangulacao, user_resolution=user_resolution)

        # Adiciona nós "fantasmas" para o D3.js (apenas para a paleta de cores funcionar)
        tree["nodes"].append({"id": "lado_A", "label": "Lado A", "side": 0})
        tree["nodes"].append({"id": "lado_B", "label": "Lado B", "side": 1})

        return render_template(
            "resultado.html", 
            tree_json=json.dumps(tree, ensure_ascii=False),
            table_data=table_data # <-- ADICIONADO
        )

    return render_template("index.html")


# ------------------------------------------------------------
# RUN
# ------------------------------------------------------------
if __name__ == "__main__":
    app.run(debug=True)