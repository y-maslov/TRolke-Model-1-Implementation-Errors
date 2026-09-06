
#include <iostream>
#include <fstream>
#include <set>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include <cstdlib>
#include <string>
#include <math.h>
#include <vector>
#include "TRolke.h"
#include "TRolkeDefault.h"
#include "TRolke128.h"
#include "TRolke128_num.h"
#include "TFeldmanCousins.h"
#include <quadmath.h>
#include <vector>
#include <utility>
#include "TPad.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TMultiGraph.h"


// rootcomp TRolke128_num.cxx TRolke128_num.h TRolke128.cxx TRolke128.h LL_comparsion.cpp -l quadmath -fext-numeric-literals -o LL_comparsion.exe

template <typename TRolkeX, typename float_type> class LL_Grapher {
    public:

        LL_Grapher(float_type alpha) {
            SetCL(alpha);
        }

        void SetCL(float_type alpha) {
            m_calculator.SetCL(alpha);
        }

        void CalcUnfixedMu(int x, int y, int z, float_type tau, int m, std::vector <double> mu_vec, float_type bias = 0) {
            m_mu_vec = mu_vec;
            m_LL_vec.resize(m_mu_vec.size(), 0.0Q);
            for (size_t i = 0; i < m_mu_vec.size(); i ++) {
                float_type mu = mu_vec[i];
                m_LL_vec[i] = - m_calculator.EvalLikeMod1(mu, x, y, z, tau, m, 3) + bias;
            }
        }

        void AddToMultiGraphFixedMu(TMultiGraph* multiGraph, TLegend* legend = nullptr, int col_bias = 0, std::string title = "title") {
            auto graph = new TGraph(m_mu_vec.size(), &m_mu_vec[0], &m_LL_vec[0]);
            graph->SetLineColor(kBlack + col_bias);
            graph->SetTitle(title.c_str());
            multiGraph->Add(graph);
            if (legend != nullptr) {
                legend->AddEntry(graph, graph->GetTitle(),"l");
            }
        }

    private:

        TRolkeX m_calculator;

        std::vector <double> m_mu_vec;
        std::vector <double> m_LL_vec;
};

void DrawLLGraphs(int x, int y, int z, int m, __float128 tau, std::vector <double> mu_vec) {
    double alpha = 0.95;

    auto legend = new TLegend(0.1, 0.9, 0.8, 0.8);

    auto mg = new TMultiGraph();
    std::string title = "x = " + std::to_string(x) + " y = "  + std::to_string(y) + 
        " z = "  + std::to_string(z) + " m = "  + std::to_string(m) + " tau = "  + std::to_string((double) tau);
    mg->SetTitle(title.c_str());

    LL_Grapher <TRolkeDefault, Double_t> defaultRolke(alpha); 
    defaultRolke.CalcUnfixedMu(x, y, z, tau, m, mu_vec);
    defaultRolke.AddToMultiGraphFixedMu(mg, legend, 0, "defaultRolke");

    // LL_Grapher <TRolke128, __float128> defaultRolke128(alpha); 
    // defaultRolke128.CalcUnfixedMu(x, y, z, tau, m, mu_vec);
    // defaultRolke128.AddToMultiGraphFixedMu(mg, legend, 1, "defaultRolke128");

    LL_Grapher <TRolke128_num, __float128> updatedRolke(alpha); 
    updatedRolke.CalcUnfixedMu(x, y, z, tau, m, mu_vec);
    updatedRolke.AddToMultiGraphFixedMu(mg, legend, 1, "updatedRolke");

    gPad->SetLogx(1);
    gPad->SetLogy(1);

    mg->Draw("APL");

    mg->GetXaxis()->SetTitle("mu");
    mg->GetYaxis()->SetTitle("- sup(2logLL(mu0, b, e); b, e)");

    gPad->Modified();

    legend->Draw();

}


int main(int argc, char **argv) {

    TApplication app("app", &argc, argv);


    std::vector <int> x_vec =           {  32,   32,      83,      83,      43,      43,         8,       3,  6700,  83000};
    std::vector <int> y_vec =           {  32,   32,      32,      32,      32,      32,         4,       3,    67,  67000};
    std::vector <int> z_vec =           { 259,  259211,  259,  259211,     259,  259211,        85,  259211,   285,    259};
    std::vector <int> m_vec =           {1000, 1000000, 1000, 1000000,    1000, 1000000,       100, 1000000,  1000,   1000};
    std::vector <__float128> tau_vec =  {   1,       1,    1,       1,       1,       1,       3.5,       1,  0.01,      1};
 
    std::vector <double> mu_vec;

    for (size_t i = 0; i < 10; i +=1) {
        mu_vec.push_back(i);
    }
    for (size_t i = 10; i < 100; i +=5) {
        mu_vec.push_back(i);
    }
    for (size_t i = 100; i < 1000; i +=50) {
        mu_vec.push_back(i);
    }
    for (size_t i = 1000; i < 10000; i +=500) {
        mu_vec.push_back(i);
    }
    for (size_t i = 10000; i < 100000; i +=5000) {
        mu_vec.push_back(i);
    }
    for (size_t i = 100000; i < 1000000; i +=50000) {
        mu_vec.push_back(i);
    }



    auto canvas = new TCanvas("canvas", "canvas", 1500, 1000);
    canvas->Divide(2);

    for (size_t i = 8; i < 10; i ++) {
        canvas->cd(i + 1 - 8);
        DrawLLGraphs(x_vec[i], y_vec[i], z_vec[i], m_vec[i], tau_vec[i], mu_vec);
    }

    app.Run();
    return 0;
}
