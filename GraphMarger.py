import ROOT

def draw_multi_graphs():
    filenames = [
        "/Users/shohtatakami/physics/COMETDS/DS189A_944turns/Diagonal_Scan_Simulation/x_neg760__x1_10.000000__x2_-10.000000_out.root",
        "/Users/shohtatakami/physics/COMETDS/DS189A_944turns/Diagonal_Scan_Simulation/x_neg760__x1_0.000000__x2_0.000000_out.root",
        "/Users/shohtatakami/physics/COMETDS/DS189A_944turns_Integration/Diagonal_Scan_Simulation/x_neg760__x1_10.000000__x2_-10.000000_out.root",
        "/Users/shohtatakami/physics/COMETDS/DS189A_944turns_Integration/Diagonal_Scan_Simulation/x_neg760__x1_0.000000__x2_0.000000_out.root"
    ]

    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta]
    
    labels = [
        "x1=10, x2=-10 (Before)",
        "x1=0, x2=0 (Before)",
        "x1=10, x2=-10 (After)",
        "x1=0, x2=0 (After)"
    ]

    # まとめて描画するためのキャンバス
    combined_canvas = ROOT.TCanvas("combined_canvas", "All Graphs Overlaid", 800, 600)
    legend = ROOT.TLegend(0.7, 0.1, 0.9, 0.3)

    # それぞれのキャンバスを保存するリスト
    individual_canvases = []

    graphs = []
    for i, filename in enumerate(filenames):
        file = ROOT.TFile.Open(filename, "READ")
        if not file or file.IsZombie():
            print(f"Error: Cannot open file {filename}")
            continue

        temp_canvas = file.Get("cxneg760")
        if not temp_canvas:
            print(f"Error: Cannot find canvas in {filename}")
            file.Close()
            continue

        graph = None
        for obj in temp_canvas.GetListOfPrimitives():
            if isinstance(obj, ROOT.TGraph):
                graph = obj
                break

        if not graph:
            print(f"Error: Cannot find TGraph in {filename}")
            file.Close()
            continue

        cloned_graph = graph.Clone(f"graph_{i}")
        cloned_graph.SetLineColor(colors[i])
        cloned_graph.SetMarkerColor(colors[i])
        cloned_graph.SetMarkerStyle(20 + i)

        # まとめたキャンバスに追加
        draw_option = "AP" if i == 0 else "P SAME"
        combined_canvas.cd()
        cloned_graph.Draw(draw_option)
        legend.AddEntry(cloned_graph, labels[i], "lp")  # **ラベルを適用**

        # 個別キャンバスの作成
        individual_canvas = ROOT.TCanvas(f"canvas_{i}", f"Graph {i+1}", 800, 600)
        individual_canvases.append(individual_canvas)
        individual_canvas.cd()
        cloned_graph.Draw("AP")
        individual_canvas.Update()

        graphs.append(cloned_graph)
        file.Close()

    # まとめたキャンバスを保存
    combined_canvas.cd()
    legend.Draw()
    combined_canvas.SetGridx()
    combined_canvas.SetGridy()
    combined_canvas.Update()
    rootplot_dir = "/Users/shohtatakami/github/COMET_DS_MFM/"
    outputFile = ROOT.TFile(f"{rootplot_dir}combined.root", "RECREATE")
    combined_canvas.Write()
    outputFile.Close()
    print(f"Successfuly Saved : {outputFile}")
    combined_canvas.SaveAs("overlayed_graphs.png")

    # 個別のキャンバスを保存
    for i, canvas in enumerate(individual_canvases):
        canvas.SaveAs(f"graph_{i+1}.png")

    print("Graphs saved successfully.")

# スクリプトを実行
draw_multi_graphs()