import ROOT as rt


def setTDRStyle():
    tdrStyle = rt.TStyle("tdrStyle", "Style for P-TDR")

    # for the canvas:
    tdrStyle.SetCanvasBorderMode(0)
    tdrStyle.SetCanvasColor(rt.kWhite)
    tdrStyle.SetCanvasDefH(600)  # Height of canvas
    tdrStyle.SetCanvasDefW(800)  # Width of canvas
    tdrStyle.SetCanvasDefX(0)  # POsition on screen
    tdrStyle.SetCanvasDefY(0)

    tdrStyle.SetPadBorderMode(0)
    # tdrStyle.SetPadBorderSize(Width_t size = 1)
    tdrStyle.SetPadColor(rt.kWhite)
    tdrStyle.SetPadGridX(False)
    tdrStyle.SetPadGridY(False)
    tdrStyle.SetGridColor(0)
    tdrStyle.SetGridStyle(3)
    tdrStyle.SetGridWidth(1)

    # For the frame:
    tdrStyle.SetFrameBorderMode(0)
    tdrStyle.SetFrameBorderSize(1)
    tdrStyle.SetFrameFillColor(0)
    tdrStyle.SetFrameFillStyle(0)
    tdrStyle.SetFrameLineColor(1)
    tdrStyle.SetFrameLineStyle(1)
    tdrStyle.SetFrameLineWidth(1)

    # For the histo:
    # tdrStyle.SetHistFillColor(1)
    # tdrStyle.SetHistFillStyle(0)
    tdrStyle.SetHistLineColor(1)
    tdrStyle.SetHistLineStyle(0)
    tdrStyle.SetHistLineWidth(1)
    # tdrStyle.SetLegoInnerR(Float_t rad = 0.5)
    # tdrStyle.SetNumberContours(Int_t number = 20)

    tdrStyle.SetEndErrorSize(2)
    # tdrStyle.SetErrorMarker(20)
    # tdrStyle.SetErrorX(0.)

    tdrStyle.SetMarkerStyle(20)

    # For the fit/function:
    tdrStyle.SetOptFit(1)
    tdrStyle.SetFitFormat("5.4g")
    tdrStyle.SetFuncColor(2)
    tdrStyle.SetFuncStyle(1)
    tdrStyle.SetFuncWidth(1)

    # For the date:
    tdrStyle.SetOptDate(0)
    # tdrStyle.SetDateX(Float_t x = 0.01)
    # tdrStyle.SetDateY(Float_t y = 0.01)

    # For the statistics box:
    tdrStyle.SetOptFile(0)
    tdrStyle.SetOptStat(0)  # To display the mean and RMS:   SetOptStat("mr")
    tdrStyle.SetStatColor(rt.kWhite)
    tdrStyle.SetStatFont(42)
    tdrStyle.SetStatFontSize(0.025)
    tdrStyle.SetStatTextColor(1)
    tdrStyle.SetStatFormat("6.4g")
    tdrStyle.SetStatBorderSize(1)
    tdrStyle.SetStatH(0.1)
    tdrStyle.SetStatW(0.15)
    # tdrStyle.SetStatStyle(Style_t style = 1001)
    # tdrStyle.SetStatX(Float_t x = 0)
    # tdrStyle.SetStatY(Float_t y = 0)

    # Margins:
    tdrStyle.SetPadTopMargin(0.05)
    tdrStyle.SetPadBottomMargin(0.13)
    tdrStyle.SetPadLeftMargin(0.16)
    tdrStyle.SetPadRightMargin(0.02)

    # For the Global title:

    tdrStyle.SetOptTitle(0)
    tdrStyle.SetTitleFont(42)
    tdrStyle.SetTitleColor(1)
    tdrStyle.SetTitleTextColor(1)
    tdrStyle.SetTitleFillColor(10)
    tdrStyle.SetTitleFontSize(0.05)
    # tdrStyle.SetTitleH(0) # Set the height of the title box
    # tdrStyle.SetTitleW(0) # Set the width of the title box
    # tdrStyle.SetTitleX(0) # Set the position of the title box
    # tdrStyle.SetTitleY(0.985) # Set the position of the title box
    # tdrStyle.SetTitleStyle(Style_t style = 1001)
    # tdrStyle.SetTitleBorderSize(2)

    # For the axis titles:

    tdrStyle.SetTitleColor(1, "XYZ")
    tdrStyle.SetTitleFont(42, "XYZ")
    tdrStyle.SetTitleSize(0.04, "XYZ")
    # tdrStyle.SetTitleXSize(Float_t size = 0.02) # Another way to set the size?
    # tdrStyle.SetTitleYSize(Float_t size = 0.02)
    tdrStyle.SetTitleXOffset(0.9)
    tdrStyle.SetTitleYOffset(1.25)
    # tdrStyle.SetTitleOffset(1.1, "Y") # Another way to set the Offset

    # For the axis labels:

    tdrStyle.SetLabelColor(1, "XYZ")
    tdrStyle.SetLabelFont(42, "XYZ")
    tdrStyle.SetLabelOffset(0.02, "XYZ")
    tdrStyle.SetLabelSize(0.04, "XYZ")

    # For the axis:

    tdrStyle.SetAxisColor(1, "XYZ")
    tdrStyle.SetStripDecimals(True)
    tdrStyle.SetTickLength(0.03, "XYZ")
    # 3    tdrStyle.SetNdivisions(510, "XYZ")
    tdrStyle.SetNdivisions(510, "XYZ")
    tdrStyle.SetNdivisions(505, "X")
    tdrStyle.SetPadTickX(1)  # To get tick marks on the opposite side of the frame
    tdrStyle.SetPadTickY(1)

    # Change for log plots:
    tdrStyle.SetOptLogx(0)
    tdrStyle.SetOptLogy(0)
    tdrStyle.SetOptLogz(0)

    # Postscript options:
    tdrStyle.SetPaperSize(20.0, 20.0)
    # tdrStyle.SetLineScalePS(Float_t scale = 3)
    # tdrStyle.SetLineStyleString(Int_t i, const char* text)
    # tdrStyle.SetHeaderPS(const char* header)
    # tdrStyle.SetTitlePS(const char* pstitle)

    # tdrStyle.SetBarOffset(Float_t baroff = 0.5)
    # tdrStyle.SetBarWidth(Float_t barwidth = 0.5)
    # tdrStyle.SetPaintTextFormat(const char* format = "g")
    # tdrStyle.SetPalette(Int_t ncolors = 0, Int_t* colors = 0)
    # tdrStyle.SetTimeOffset(Double_t toffset)
    # tdrStyle.SetHistMinimumZero(kTRUE)

    tdrStyle.SetHatchesLineWidth(5)
    tdrStyle.SetHatchesSpacing(0.05)

    tdrStyle.cd()


# CMS_lumi
#   Initiated by: Gautier Hamel de Monchenault (Saclay)
#   Translated in Python by: Joshua Hardenbrook (Princeton)
#   Updated by:   Dinko Ferencek (Rutgers)
#

cmsText = "CMS"
cmsTextFont = 61
cmsTextConsumedWidth = 0.10

writeExtraText = True
# extraText = "Simulation Preliminary"
extraText = "Preliminary"
extraTextFont = 52

lumiTextSize = 0.6
lumiTextOffset = 0.2

cmsTextSize = 0.8
cmsTextOffset = 0.1

relPosX = 0.045
relPosY = 0.035
relExtraDY = 1.2

extraOverCmsTextSize = 0.76

lumi_13TeV = "20.1 fb^{-1}"
lumi_8TeV = "19.6 fb^{-1}"
lumi_7TeV = "5.1 fb^{-1}"
lumi_sqrtS = "138 fb^{-1} (13 TeV)"

drawLogo = False


def CMS_lumi(pad, iPeriod, iPosX):
    outOfFrame = False
    if iPosX / 10 == 0:
        outOfFrame = True

    alignY_ = 3
    alignX_ = 2
    if iPosX / 10 == 0:
        alignX_ = 1
    if iPosX == 0:
        alignY_ = 1
    if iPosX / 10 == 1:
        alignX_ = 1
    if iPosX / 10 == 2:
        alignX_ = 2
    if iPosX / 10 == 3:
        alignX_ = 3
    align_ = 10 * alignX_ + alignY_

    H = pad.GetWh()
    W = pad.GetWw()
    l = pad.GetLeftMargin()
    t = pad.GetTopMargin()
    r = pad.GetRightMargin() + 0.05
    b = pad.GetBottomMargin()

    pad.cd()

    lumiText = ""
    if iPeriod == 1:
        lumiText += lumi_7TeV
        lumiText += " (7 TeV)"
    elif iPeriod == 2:
        lumiText += lumi_8TeV
        lumiText += " (8 TeV)"

    elif iPeriod == 3:
        lumiText = lumi_8TeV
        lumiText += " (8 TeV)"
        lumiText += " + "
        lumiText += lumi_7TeV
        lumiText += " (7 TeV)"
    elif iPeriod == 4:
        lumiText += lumi_13TeV
        lumiText += " (13 TeV)"
    elif iPeriod == 7:
        if outOfFrame:
            lumiText += "#scale[0.85]{"
        lumiText += lumi_13TeV
        lumiText += " (13 TeV)"
        lumiText += " + "
        lumiText += lumi_8TeV
        lumiText += " (8 TeV)"
        lumiText += " + "
        lumiText += lumi_7TeV
        lumiText += " (7 TeV)"
        if outOfFrame:
            lumiText += "}"
    elif iPeriod == 12:
        lumiText += "8 TeV"
    elif iPeriod == 0:
        lumiText += lumi_sqrtS

    print(lumiText)

    latex = rt.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(rt.kBlack)

    extraTextSize = extraOverCmsTextSize * cmsTextSize

    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    latex.SetTextSize(lumiTextSize * t)
    print(1 - r, 1 - t + lumiTextOffset * t, lumiText)
    latex.DrawLatex(1 - r - 0.03, 1 - t + lumiTextOffset * t, lumiText)

    if outOfFrame:
        latex.SetTextFont(cmsTextFont)
        latex.SetTextAlign(11)
        latex.SetTextSize(cmsTextSize * t)
        latex.DrawLatex(l, 1 - t + lumiTextOffset * t, cmsText)

    pad.cd()

    posX_ = 0
    if iPosX % 10 <= 1:
        posX_ = l + relPosX * (1 - l - r)
    elif iPosX % 10 == 2:
        posX_ = l + 0.5 * (1 - l - r)
    elif iPosX % 10 == 3:
        posX_ = 1 - r - relPosX * (1 - l - r)

    posY_ = 1 - t - relPosY * (1 - t - b)
    print("posX_", posX_)
    posX_ = 0.05
    if not outOfFrame:
        if drawLogo:
            posX_ = l + 0.045 * (1 - l - r) * W / H
            posY_ = 1 - t - 0.045 * (1 - t - b)
            xl_0 = posX_
            yl_0 = posY_ - 0.15
            xl_1 = posX_ + 0.15 * H / W
            yl_1 = posY_
            CMS_logo = rt.TASImage("CMS-BW-label.png")
            pad_logo = rt.TPad("logo", "logo", xl_0, yl_0, xl_1, yl_1)
            pad_logo.Draw()
            pad_logo.cd()
            CMS_logo.Draw("X")
            pad_logo.Modified()
            pad.cd()
        else:
            latex.SetTextFont(cmsTextFont)
            latex.SetTextSize(cmsTextSize * t)
            latex.SetTextAlign(align_)
            posX_ = 0.1
            latex.DrawLatex(posX_, posY_, cmsText)
            if writeExtraText:
                latex.SetTextFont(extraTextFont)
                latex.SetTextAlign(align_)
                latex.SetTextSize(extraTextSize * t)
                latex.DrawLatex(posX_, posY_ - relExtraDY * cmsTextSize * t, extraText)
    elif writeExtraText:
        if iPosX == 0:
            posX_ = l + relPosX * (1 - l - r) + cmsTextConsumedWidth
            posY_ = 1 - t + lumiTextOffset * t
        latex.SetTextFont(extraTextFont)
        latex.SetTextSize(extraTextSize * t)
        latex.SetTextAlign(align_)
        latex.DrawLatex(posX_, posY_, extraText)

    pad.Update()
