// ======================================================
//       D3 — Árvore Vertical Avançada ArboreDNA (corrigido)
// ======================================================

(function () {

    const container = d3.select("#tree-container");
    const width = container.node().clientWidth;
    const height = container.node().clientHeight;

    // --- MUDANÇA 1: Criar o handler de zoom ---
    const zoomHandler = d3.zoom()
        .scaleExtent([0.3, 2.5])
        .on("zoom", (event) => g.attr("transform", event.transform));

    const svg = container.append("svg")
        .attr("width", width)
        .attr("height", height)
        .style("cursor", "grab")
        .call(zoomHandler); // --- MUDANÇA 2: Chamar o handler

    const g = svg.append("g").attr("transform", "translate(60,60)");

    // safety checks
    if (typeof DATA === "undefined" || !DATA || !Array.isArray(DATA.nodes) || !Array.isArray(DATA.edges)) {
        // show a message inside the container
        container.append("div")
            .style("padding", "20px")
            .style("color", "#b00020")
            .text("Dados inválidos para a árvore (DATA não encontrado ou formato incorreto).");
        return;
    }

    // ---------------------------------------------
    // 1. Converter nodes + edges → hierarquia D3
    // ---------------------------------------------
    const nodesById = {};
    DATA.nodes.forEach(n => nodesById[n.id] = { ...n, children: [], _children: [] });

    DATA.edges.forEach(e => {
        const parent = nodesById[e.from];
        const child  = nodesById[e.to];
        if (parent && child) parent.children.push(child);
    });

    if (!nodesById["root"]) {
        container.append("div")
            .style("padding", "20px")
            .style("color", "#b00020")
            .text("Raiz 'root' não encontrada nos nodes.");
        return;
    }

    const root = d3.hierarchy(nodesById["root"]);

    // ---------------------------------------------
    // 2. Configurar layout da árvore
    // ---------------------------------------------
    const treeLayout = d3
        .tree()
        .nodeSize([90, 140]); // espaçamento horizontal (x) e vertical (y)

    root.x0 = width / 2;
    root.y0 = 0;

    // **Observação**: não colapsar tudo por padrão — queremos ver os matches.
    // Se preferir colapsar, descomente a linha abaixo:
    // root.children?.forEach(collapse);

    update(root);

    // ---------------------------------------------
    // Collapse recursivo
    // ---------------------------------------------
    function collapse(node) {
        if (node.children) {
            node._children = node.children;
            node._children.forEach(collapse);
            node.children = null;
        }
    }

    // ---------------------------------------------
    // Função de update principal
    // ---------------------------------------------
    function update(source) {

		// +++ NOVA ESCALA DE TAMANHO +++
		// Encontra o max cM (em nós folha) para a escala
		const maxCm = root.descendants().reduce((max, d) => {
			// Apenas matches (que têm 'cm' e não têm 'children')
			if (d.data.cm && !d.children) {
				return Math.max(max, d.data.cm);
			}
			return max;
		}, 70); // min

		// Escala: de 70cM (min) até o maxCm
		// Mapeia para um raio de 10px (pequeno) até 20px (grande)
		const radiusScale = d3.scaleSqrt()
			.domain([70, maxCm])
			.range([10, 20])
			.clamp(true); // não deixa passar de 20 ou ser menor que 10
		// +++ FIM DA NOVA ESCALA +++

		treeLayout(root);

        // --------------------------
        // LINKS
        // --------------------------
        const link = g.selectAll("path.link")
            .data(root.links(), d => d.target.data.id);

        // ENTER
        const linkEnter = link.enter()
            .append("path")
            .attr("class", "link")
            .attr("fill", "none")
            .attr("stroke", "#b6b6b6")
            .attr("stroke-width", 2)
            .attr("d", d => diagonal({ source: source, target: source }));

        // UPDATE + MERGE
        linkEnter.merge(link)
            .transition().duration(450)
            .attr("d", diagonal);

        // EXIT
        link.exit()
            .transition().duration(300)
            .attr("d", d => diagonal({ source: source, target: source }))
            .remove();

        // --------------------------
        // NODES
        // --------------------------
        const node = g.selectAll("g.node")
            .data(root.descendants(), d => d.data.id);

        const nodeEnter = node.enter()
            .append("g")
            .attr("class", "node")
            .attr("transform", d => `translate(${source.x0}, ${source.y0})`)
            .style("cursor", "pointer")
            // D3 v7: event is first param, datum second
            .on("click", (event, d) => toggle(d))
            .on("mouseover", (event, d) => showTooltip(event, d))
            .on("mouseout", hideTooltip);

        nodeEnter.append("circle")
            .attr("r", 1e-6)
            .attr("fill", d => color(d.data.side))
            .attr("stroke", "#333")
            .attr("stroke-width", 1.5);

        // --- CORREÇÃO DE LEITURA: ROTAÇÃO DE TEXTO ---
        nodeEnter.append("text")
            .attr("font-family", "sans-serif")
            .attr("font-size", "11px")
            .attr("fill", "#222")
            .text(d => d.data.label)
            
            // LÓGICA DE POSIÇÃO E ROTAÇÃO
            .attr("transform", d => {
                // SE FOR MATCH (Folha/Sem filhos):
                if (!d.children && !d._children) {
                    // Move 20px para baixo e rotaciona 90 graus
                    return "translate(0, 20) rotate(90)";
                }
                // SE FOR GRUPO (Tem filhos):
                // Mantém horizontal e move 20px para cima
                return "translate(0, -20)";
            })
            
            // LÓGICA DE ALINHAMENTO
            .attr("text-anchor", d => {
                // Match: "start" faz o texto começar na bolinha e descer
                if (!d.children && !d._children) return "start"; 
                // Grupo: "middle" centraliza o texto sobre a bolinha
                return "middle";
            })

            .style("opacity", 0)
            // Sombra branca para leitura
            .style("text-shadow", "2px 2px 4px white, -2px -2px 4px white, -2px 2px 4px white, 2px -2px 4px white");

        // UPDATE
        const nodeUpdate = nodeEnter.merge(node);

        nodeUpdate.transition().duration(450)
            .attr("transform", d => `translate(${d.x}, ${d.y})`);

        nodeUpdate.select("circle")
            .transition().duration(450)
			// +++ ATRIBUTO 'r' ATUALIZADO +++
			.attr("r", d => {
				if (d.data.id === "root") return 18; // Raiz grande
				if (d.data.cm && !d.children) {
					return radiusScale(d.data.cm); // Tamanho baseado em cM
				}
				return 16; // Tamanho padrão (Subramos, etc)
			})
			// +++ FIM DA ATUALIZAÇÃO +C
            .attr("fill", d => color(d.data.side));

        nodeUpdate.select("text")
            .transition().duration(450)
            .style("opacity", 1);

        // EXIT
        const nodeExit = node.exit()
            .transition().duration(300)
            .attr("transform", d => `translate(${source.x}, ${source.y})`)
            .remove();

        nodeExit.select("circle").attr("r", 1e-6);
        nodeExit.select("text").style("opacity", 0);

        // Atualizar antigas posições
        root.descendants().forEach(d => {
            d.x0 = d.x;
            d.y0 = d.y;
        });
    }

    // ---------------------------------------------
    // Toggle de expandir/colapsar
    // ---------------------------------------------
    function toggle(d) {
        if (d.children) {
            d._children = d.children;
            d.children = null;
        } else {
            d.children = d._children;
            d._children = null;
        }
        update(d);
    }

    // ---------------------------------------------
    // Tooltip (D3 v7: receber event)
    // ---------------------------------------------
    const tooltip = d3.select("body")
        .append("div")
        .style("position", "absolute")
        .style("padding", "6px 12px")
        .style("background", "#ffffffdd")
        .style("border-radius", "8px")
        .style("border", "1px solid #777")
        .style("font-size", "13px")
        .style("pointer-events", "none")
        .style("opacity", 0);

    function showTooltip(event, d) {
        if (!d || !d.data) return;
        const label = d.data.label || "";
        const sideText = d.data.side === 0 ? "Lado A" : (d.data.side === 1 ? "Lado B" : "");

        tooltip
            .style("opacity", 1)
            .html(`<strong>${label}</strong><br>${sideText}`)
            .style("left", (event.pageX + 15) + "px")
            .style("top", (event.pageY + 15) + "px");
    }

    function hideTooltip() {
        tooltip.style("opacity", 0);
    }

    // ---------------------------------------------
    // Linha curva vertical
    // ---------------------------------------------
    function diagonal(d) {
        return `
            M ${d.source.x},${d.source.y}
            C ${d.source.x}, ${(d.source.y + d.target.y) / 2}
              ${d.target.x}, ${(d.source.y + d.target.y) / 2}
              ${d.target.x}, ${d.target.y}
        `;
    }

    // ---------------------------------------------
    // Cor moderna (MODO ENDOGAMIA)
    // ---------------------------------------------
    function color(side_index) {
        // Mapeia o 'side' (que agora é um índice) para cores
        const colors = [
            "#F44336", // 0 = A (vermelho)
            "#2196F3", // 1 = B (azul)
            "#4CAF50", // 2 (verde)
            "#FFC107", // 3 (amarelo)
            "#9C27B0", // 4 (roxo)
            "#00BCD4", // 5 (ciano)
            "#E91E63"  // 6 (rosa)
        ];
        
        if (typeof side_index === "number" && side_index >= 0) {
            return colors[side_index % colors.length];
        }
        return "#9E9E9E"; // Cor Padrão (cinza)
    }

    // ---------------------------------------------
    // --- MUDANÇA 3: Adicionar listeners aos botões ---
    // ---------------------------------------------
    d3.select("#zoom-in").on("click", function() {
        svg.transition().duration(250).call(zoomHandler.scaleBy, 1.3); // Zoom in 30%
    });

    d3.select("#zoom-out").on("click", function() {
        svg.transition().duration(250).call(zoomHandler.scaleBy, 0.7); // Zoom out 30%
    });
	
})();