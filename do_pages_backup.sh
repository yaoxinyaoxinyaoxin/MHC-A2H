#!/bin/bash
cd "/Users/yaoxin/Library/Mobile Documents/com~apple~CloudDocs/AI_Code_yun/project/8_多组学与带状疱疹/0_汇总所有脚本(代码可用性_github)/4_进行github化/3_github化(英文注释版本)/Github_Scripts_English_Only_20260504_190907（来自备份）"

git checkout --orphan gh-pages
git rm -rf .
cp -r "/Users/yaoxin/Library/Mobile Documents/com~apple~CloudDocs/AI_Code_yun/project/8_多组学与带状疱疹/0_汇总所有脚本(代码可用性_github)/4_进行github化/3_github化(英文注释版本)/3D/3D_edge-bundling_plot_Interactive" ./
cp -r "/Users/yaoxin/Library/Mobile Documents/com~apple~CloudDocs/AI_Code_yun/project/8_多组学与带状疱疹/0_汇总所有脚本(代码可用性_github)/4_进行github化/3_github化(英文注释版本)/3D/3D_Scatter_Plot_Interactive" ./

cat << 'EOF' > index.html
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MHC-A2H Interactive 3D Plots</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }
        h1 { color: #333; }
        ul { list-style-type: none; padding: 0; }
        li { margin: 10px 0; }
        a { text-decoration: none; color: #0066cc; font-size: 18px; }
        a:hover { text-decoration: underline; }
    </style>
</head>
<body>
    <h1>MHC-A2H Interactive 3D Plots</h1>
    <p>Welcome to the interactive data visualization page for our manuscript. Below are the 3D plots available for exploration:</p>
    <ul>
        <li><a href="3D_edge-bundling_plot_Interactive/3D_Interactive_Network.html" target="_blank">🔍 Supplementary File 2: 3D Interactive Causal Network (Hierarchical Edge-Bundling)</a></li>
        <li><a href="3D_Scatter_Plot_Interactive/3D_Scatter_Plot_Interactive_Complete.html" target="_blank">🔍 Supplementary File 1: 3D Scatter Plot</a></li>
    </ul>
    <p><br><em>Note: These plots are fully interactive. You can rotate, zoom, and hover over data points for more details.</em></p>
</body>
</html>
EOF

git add .
git commit -m "Add 3D interactive plots for GitHub Pages"
git push -f -u origin gh-pages
git checkout main
