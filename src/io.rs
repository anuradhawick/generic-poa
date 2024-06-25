use crate::graph::{EdgeData, NodeData};
use petgraph::{
    dot::{Config, Dot},
    graph::DiGraph,
    Graph,
};
use std::{
    cmp::min,
    fs::File,
    io::{BufWriter, Write},
};

fn get_dot(graph: &DiGraph<NodeData, EdgeData>) -> Dot<&Graph<NodeData, EdgeData>> {
    Dot::with_attr_getters(
        graph,
        &[Config::EdgeNoLabel, Config::NodeNoLabel],
        &|_, e| {
            format!(
                "label=\"Fragments: [{}]\" penwidth={1} minlen={1}",
                e.weight().labels.join(", "),
                min(10, e.weight().labels.len())
            )
        },
        &|_, e| format!("label = \"{}\"", e.1.item),
    )
}

pub fn write_dot(graph: &DiGraph<NodeData, EdgeData>, path: &str) -> Result<(), String> {
    let dot = get_dot(graph);
    let file = File::create(format!("{}.graph.dot", path))
        .map_err(|_| format!("Unable to create file: {}", path))?;
    let mut writer = BufWriter::new(file);
    writer
        .write_all(format!("{:?}", dot).as_bytes())
        .map_err(|_| "IO Error".to_string())
}

pub fn write_html(graph: &DiGraph<NodeData, EdgeData>, path: &str) -> Result<(), String> {
    let dot = get_dot(graph);
    let html = format!(
        r#"<!DOCTYPE html>
<html lang="en">
<head>
    <title>Generic POA Graph</title>
    <script src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"></script>
    <style>
        #mynetwork {{
            height: 100vh;
        }}
    </style>
</head>
<body>
    <div id="mynetwork" style></div>

    <script type="text/javascript">
    var container = document.getElementById("mynetwork");
    var dot = `{:?}`;
    var data = vis.parseDOTNetwork(dot);
    var network = new vis.Network(container, data);
    </script>
</body>
</html>
"#,
        dot
    );
    let file = File::create(format!("{}.graph.html", path))
        .map_err(|_| format!("Unable to create file: {}", path))?;
    let mut writer = BufWriter::new(file);
    writer
        .write_all(html.as_bytes())
        .map_err(|_| "IO Error".to_string())
}
