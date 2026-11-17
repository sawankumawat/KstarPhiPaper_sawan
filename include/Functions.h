#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "TF1.h"
#include "TMath.h"
#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include "TMinuit.h"

// ============================================================
// Power Law
// ============================================================

TF1 *PowerLaw(const Char_t *name, Double_t norm = 1.e-6,
              Double_t exponent = 6.0)
{
    TF1 *f = new TF1(name, "[0]*TMath::Power(x/5.,-[1])", 0.1, 100.);
    f->SetParameters(norm, exponent);
    f->SetParNames("Norm.", "Exponent");
    return f;
}

// ============================================================
// Levy-Tsallis
// ============================================================

Double_t FuncLevy_Sawan(Double_t *x, Double_t *par)
{
    Double_t p = (par[0] - 1) * (par[0] - 2) * par[1] * x[0] / (((pow((1 + (((sqrt((par[2] * par[2]) + (x[0] * x[0]))) - par[2]) / (par[0] * par[3]))), par[0]) * (par[0] * par[3] * ((par[0] * par[3]) + (par[2] * (par[0] - 2)))))));
    return (p);
}

Double_t LevyTsallis_Func(const Double_t *x, const Double_t *p)
{
    Double_t pt = x[0], mass = p[0], mt = TMath::Sqrt(pt * pt + mass * mass), n = p[1], C = p[2], norm = p[3];
    Double_t part1 = (n - 1.) * (n - 2.), part2 = n * C * (n * C + mass * (n - 2.)), part3 = part1 / part2, part4 = 1. + (mt - mass) / (n * C);
    return pt * norm * part3 * TMath::Power(part4, -n);
}

TF1 *LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.)
{
    TF1 *f = new TF1(name, LevyTsallis_Func, 0., 25., 4);
    f->SetParameters(mass, n, C, norm);
    f->SetParNames("mass", "n", "C", "norm");
    f->FixParameter(0, mass);
    f->SetParLimits(1, 2., 30.);
    f->SetParLimits(2, 0.05, 1.);
    f->SetParLimits(3, 1.e-6, 1.e6);
    return f;
}

// ============================================================
// Boltzmann
// ============================================================

Double_t Boltzmann_Func(const Double_t *x, const Double_t *p)
{
    Double_t pt = x[0], mass = p[0], mt = TMath::Sqrt(pt * pt + mass * mass), T = p[1], norm = p[2];
    return pt * norm * mt * TMath::Exp(-mt / T);
}

TF1 *Boltzmann(const Char_t *name, Double_t mass, Double_t T = 0.1, Double_t norm = 1.)
{
    TF1 *f = new TF1(name, Boltzmann_Func, 0., 20., 3);
    f->SetParameters(mass, T, norm);
    f->SetParNames("mass", "T", "norm");
    f->FixParameter(0, mass);
    return f;
}

/*****************************************************************/
/* BOLTZMANN-GIBBS BLAST-WAVE */
/*****************************************************************/

static TF1 *fBGBlastWave_Integrand = NULL;
static TF1 *fBGBlastWave_Integrand_num = NULL;
static TF1 *fBGBlastWave_Integrand_den = NULL;
Double_t
BGBlastWave_Integrand(const Double_t *x, const Double_t *p)
{

    /*
       x[0] -> r (radius)
       p[0] -> mT (transverse mass)
       p[1] -> pT (transverse momentum)
       p[2] -> beta_max (surface velocity)
       p[3] -> T (freezout temperature)
       p[4] -> n (velocity profile)
    */

    Double_t r = x[0];
    Double_t mt = p[0];
    Double_t pt = p[1];
    Double_t beta_max = p[2];
    Double_t temp_1 = 1. / p[3];
    Double_t n = p[4];

    Double_t beta = beta_max * TMath::Power(r, n);
    if (beta > 0.9999999999999999)
        beta = 0.9999999999999999;
    Double_t rho = TMath::ATanH(beta);
    Double_t argI0 = pt * TMath::SinH(rho) * temp_1;
    if (argI0 > 700.)
        argI0 = 700.;
    Double_t argK1 = mt * TMath::CosH(rho) * temp_1;
    //  if (argI0 > 100 || argI0 < -100)
    //    printf("r=%f, pt=%f, beta_max=%f, temp=%f, n=%f, mt=%f, beta=%f, rho=%f, argI0=%f, argK1=%f\n", r, pt, beta_max, 1. / temp_1, n, mt, beta, rho, argI0, argK1);
    return r * mt * TMath::BesselI0(argI0) * TMath::BesselK1(argK1);
}

Double_t
BGBlastWave_Func(const Double_t *x, const Double_t *p)
{
    /* dN/dpt */

    Double_t pt = x[0];
    Double_t mass = p[0];
    Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
    Double_t beta_max = p[1];
    Double_t temp = p[2];
    Double_t n = p[3];
    Double_t norm = p[4];

    if (!fBGBlastWave_Integrand)
        fBGBlastWave_Integrand = new TF1("fBGBlastWave_Integrand", BGBlastWave_Integrand, 0., 1., 5);
    fBGBlastWave_Integrand->SetParameters(mt, pt, beta_max, temp, n);
    Double_t integral = fBGBlastWave_Integrand->Integral(0., 1., 1.e-6);
    return norm * pt * integral;
}

Double_t
BGBlastWaveRatio_Func(const Double_t *x, const Double_t *p)
{
    /* dN/dpt */

    Double_t pt = x[0];
    Double_t mass = p[0];
    Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
    Double_t beta_max_num = p[1];
    Double_t temp_num = p[2];
    Double_t n_num = p[3];
    Double_t norm_num = p[4];
    Double_t beta_max_den = p[5];
    Double_t temp_den = p[6];
    Double_t n_den = p[7];
    Double_t norm_den = p[8];

    if (!fBGBlastWave_Integrand_num)
        fBGBlastWave_Integrand_num = new TF1("fBGBlastWave_Integrand_num", BGBlastWave_Integrand, 0., 1., 5);
    fBGBlastWave_Integrand_num->SetParameters(mt, pt, beta_max_num, temp_num, n_num);
    Double_t integral_num = fBGBlastWave_Integrand_num->Integral(0., 1.);

    if (!fBGBlastWave_Integrand_den)
        fBGBlastWave_Integrand_den = new TF1("fBGBlastWave_Integrand_den", BGBlastWave_Integrand, 0., 1., 5);
    fBGBlastWave_Integrand_den->SetParameters(mt, pt, beta_max_den, temp_den, n_den);
    Double_t integral_den = fBGBlastWave_Integrand_den->Integral(0., 1.);

    return (norm_num / norm_den) * (integral_num / integral_den);
}

Double_t
BGBlastWaveParticleRatio_Func(const Double_t *x, const Double_t *p)
{
    /* dN/dpt */

    Double_t pt = x[0];
    Double_t mass_num = p[0];
    Double_t mass_den = p[1];
    Double_t mt_num = TMath::Sqrt(pt * pt + mass_num * mass_num);
    Double_t mt_den = TMath::Sqrt(pt * pt + mass_den * mass_den);
    Double_t beta_max = p[2];
    Double_t temp = p[3];
    Double_t n = p[4];
    Double_t norm_num = p[5];
    Double_t norm_den = p[6];

    if (!fBGBlastWave_Integrand_num)
        fBGBlastWave_Integrand_num = new TF1("fBGBlastWave_Integrand_num", BGBlastWave_Integrand, 0., 1., 5);
    fBGBlastWave_Integrand_num->SetParameters(mt_num, pt, beta_max, temp, n);
    Double_t integral_num = fBGBlastWave_Integrand_num->Integral(0., 1.);

    if (!fBGBlastWave_Integrand_den)
        fBGBlastWave_Integrand_den = new TF1("fBGBlastWave_Integrand_den", BGBlastWave_Integrand, 0., 1., 5);
    fBGBlastWave_Integrand_den->SetParameters(mt_den, pt, beta_max, temp, n);
    Double_t integral_den = fBGBlastWave_Integrand_den->Integral(0., 1.);

    return (norm_num / norm_den) * (integral_num / integral_den);
}

Double_t
BGBlastWave_Func_OneOverPt(const Double_t *x, const Double_t *p)
{
    /* 1/pt dN/dpt */

    Double_t pt = x[0];
    Double_t mass = p[0];
    Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
    Double_t beta_max = p[1];
    Double_t temp = p[2];
    Double_t n = p[3];
    Double_t norm = p[4];

    if (!fBGBlastWave_Integrand)
        fBGBlastWave_Integrand = new TF1("fBGBlastWave_Integrand", BGBlastWave_Integrand, 0., 1., 5);
    fBGBlastWave_Integrand->SetParameters(mt, pt, beta_max, temp, n);
    Double_t integral = fBGBlastWave_Integrand->Integral(0., 1., 1.e-3);

    return norm * integral;
}

TF1 *BGBlastWave(const Char_t *name, Double_t mass, Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t norm = 1.e6)
{

    TF1 *fBGBlastWave = new TF1(name, BGBlastWave_Func, 0., 10., 5);
    fBGBlastWave->SetParameters(mass, beta_max, temp, n, norm);
    fBGBlastWave->SetParNames("mass", "beta_max", "T", "n", "norm");
    fBGBlastWave->FixParameter(0, mass);
    fBGBlastWave->SetParLimits(1, 0.01, 0.99);
    fBGBlastWave->SetParLimits(2, 0.01, 1.);
    fBGBlastWave->SetParLimits(3, 0.01, 50.);
    return fBGBlastWave;
}

TF1 *BGBlastWaveRatio(const Char_t *name, Double_t mass, Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t norm = 1.e6)
{

    TF1 *fBGBlastWave = new TF1(name, BGBlastWaveRatio_Func, 0., 10., 9);
    fBGBlastWave->SetParameters(mass, beta_max, temp, n, norm, beta_max, temp, n, norm);
    fBGBlastWave->SetParNames("mass", "beta_max_num", "T_num", "n_num", "norm_num", "beta_max_den", "T_den", "n_den", "norm_den");
    fBGBlastWave->FixParameter(0, mass);
    fBGBlastWave->SetParLimits(1, 0.01, 0.99);
    fBGBlastWave->SetParLimits(2, 0.01, 1.);
    fBGBlastWave->SetParLimits(3, 0.01, 10.);
    fBGBlastWave->SetParLimits(5, 0.01, 0.99);
    fBGBlastWave->SetParLimits(6, 0.01, 1.);
    fBGBlastWave->SetParLimits(7, 0.01, 10.);
    return fBGBlastWave;
}

TF1 *BGBlastWaveParticleRatio(const Char_t *name, Double_t mass_num, Double_t mass_den, Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t norm_num = 1.e6, Double_t norm_den = 1.e6)
{

    TF1 *fBGBlastWave = new TF1(name, BGBlastWaveParticleRatio_Func, 0., 10., 7);
    fBGBlastWave->SetParameters(mass_num, mass_den, beta_max, temp, n, norm_num, norm_den);
    fBGBlastWave->SetParNames("mass_num", "mass_den", "beta_max", "T", "n", "norm_num", "norm_den");
    fBGBlastWave->FixParameter(0, mass_num);
    fBGBlastWave->FixParameter(1, mass_den);
    fBGBlastWave->SetParLimits(2, 0.01, 0.99);
    fBGBlastWave->SetParLimits(3, 0.01, 1.);
    fBGBlastWave->SetParLimits(4, 0.01, 10.);
    return fBGBlastWave;
}

TF1 *BGBlastWave_OneOverPT(const Char_t *name, Double_t mass, Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t norm = 1.e6)
{

    TF1 *fBGBlastWave = new TF1(name, BGBlastWave_Func_OneOverPt, 0., 10., 5);
    fBGBlastWave->SetParameters(mass, beta_max, temp, n, norm);
    fBGBlastWave->SetParNames("mass", "beta_max", "T", "n", "norm");
    fBGBlastWave->FixParameter(0, mass);
    fBGBlastWave->SetParLimits(1, 0.01, 0.99);
    fBGBlastWave->SetParLimits(2, 0.01, 1.);
    fBGBlastWave->SetParLimits(3, 0.01, 50.);
    return fBGBlastWave;
}

// ============================================================
// mT exponential
// ============================================================

Double_t mTExp_Func(const Double_t *x, const Double_t *p)
{
    Double_t pt = x[0], mass = p[0], T = p[1], norm = p[2], mt = TMath::Sqrt(pt * pt + mass * mass);
    return pt * norm * TMath::Exp(-(mt - mass) / T);
}

TF1 *mTExp(const Char_t *name, Double_t mass, Double_t T = 0.1, Double_t norm = 1.)
{
    TF1 *f = new TF1(name, mTExp_Func, 0., 20., 3);
    f->SetParameters(mass, T, norm);
    f->SetParNames("mass", "T", "norm");
    f->FixParameter(0, mass);
    return f;
}

// ============================================================
// pT exponential
// ============================================================

Double_t pTExp_Func(const Double_t *x, const Double_t *p)
{
    Double_t pt = x[0], T = p[0], norm = p[1];
    return pt * norm * TMath::Exp(-pt / T);
}

TF1 *pTExp(const Char_t *name, Double_t T = 0.1, Double_t norm = 1.)
{
    TF1 *f = new TF1(name, pTExp_Func, 0., 20., 2);
    f->SetParameters(T, norm);
    f->SetParNames("T", "norm");
    return f;
}

// ============================================================
// Bose-Einstein
// ============================================================

Double_t BoseEinstein_Func(const Double_t *x, const Double_t *p)
{
    Double_t pt = x[0], mass = p[0], T = p[1], norm = p[2], mt = TMath::Sqrt(pt * pt + mass * mass);
    Double_t arg = mt / T;
    if (arg > 700.)
        return 0.;
    return pt * norm / (TMath::Exp(arg) - 1.);
}

TF1 *BoseEinstein(const Char_t *name, Double_t mass, Double_t T = 0.1, Double_t norm = 1.)
{
    TF1 *f = new TF1(name, BoseEinstein_Func, 0., 20., 3);
    f->SetParameters(mass, T, norm);
    f->SetParNames("mass", "T", "norm");
    f->FixParameter(0, mass);
    return f;
}

// ============================================================
// Fermi-Dirac
// ============================================================

Double_t FermiDirac_Func(const Double_t *x, const Double_t *p)
{
    Double_t pt = x[0], mass = p[0], T = p[1], norm = p[2], mt = TMath::Sqrt(pt * pt + mass * mass);
    Double_t arg = mt / T;
    if (arg > 700.)
        return 0.;
    return pt * norm / (TMath::Exp(arg) + 1.);
}

TF1 *FermiDirac(const Char_t *name, Double_t mass, Double_t T = 0.1, Double_t norm = 1.)
{
    TF1 *f = new TF1(name, FermiDirac_Func, 0., 20., 3);
    f->SetParameters(mass, T, norm);
    f->SetParNames("mass", "T", "norm");
    f->FixParameter(0, mass);
    return f;
}

/*****************************************************************/
/* TSALLIS BLAST-WAVE */
/*****************************************************************/

static TF1 *fTsallisBlastWave_Integrand_r = NULL;
Double_t
TsallisBlastWave_Integrand_r(const Double_t *x, const Double_t *p)
{
    /*
       x[0] -> r (radius)
       p[0] -> mT (transverse mass)
       p[1] -> pT (transverse momentum)
       p[2] -> beta_max (surface velocity)
       p[3] -> T (freezout temperature)
       p[4] -> n (velocity profile)
       p[5] -> q
       p[6] -> y (rapidity)
       p[7] -> phi (azimuthal angle)
    */

    Double_t r = x[0];
    Double_t mt = p[0];
    Double_t pt = p[1];
    Double_t beta_max = p[2];
    Double_t temp_1 = 1. / p[3];
    Double_t n = p[4];
    Double_t q = p[5];
    Double_t y = p[6];
    Double_t phi = p[7];

    if (q <= 1.)
        return r;

    Double_t beta = beta_max * TMath::Power(r, n);
    Double_t rho = TMath::ATanH(beta);

    Double_t part1 = mt * TMath::CosH(y) * TMath::CosH(rho);
    Double_t part2 = pt * TMath::SinH(rho) * TMath::Cos(phi);
    Double_t part3 = part1 - part2;
    Double_t part4 = 1 + (q - 1.) * temp_1 * part3;
    Double_t expo = -1. / (q - 1.);
    //  printf("part1=%f, part2=%f, part3=%f, part4=%f, expo=%f\n", part1, part2, part3, part4, expo);
    Double_t part5 = TMath::Power(part4, expo);

    return r * part5;
}

static TF1 *fTsallisBlastWave_Integrand_phi = NULL;
Double_t
TsallisBlastWave_Integrand_phi(const Double_t *x, const Double_t *p)
{
    /*
       x[0] -> phi (azimuthal angle)
    */

    Double_t phi = x[0];
    fTsallisBlastWave_Integrand_r->SetParameter(7, phi);
    Double_t integral = fTsallisBlastWave_Integrand_r->Integral(0., 1.);
    return integral;
}

static TF1 *fTsallisBlastWave_Integrand_y = NULL;
Double_t
TsallisBlastWave_Integrand_y(const Double_t *x, const Double_t *p)
{
    /*
       x[0] -> y (rapidity)
    */

    Double_t y = x[0];
    fTsallisBlastWave_Integrand_r->SetParameter(6, y);
    Double_t integral = fTsallisBlastWave_Integrand_phi->Integral(-TMath::Pi(), TMath::Pi());
    return TMath::CosH(y) * integral;
}

Double_t
TsallisBlastWave_Func(const Double_t *x, const Double_t *p)
{
    /* dN/dpt */

    Double_t pt = x[0];
    Double_t mass = p[0];
    Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
    Double_t beta_max = p[1];
    Double_t temp = p[2];
    Double_t n = p[3];
    Double_t q = p[4];
    Double_t norm = p[5];

    if (!fTsallisBlastWave_Integrand_r)
        fTsallisBlastWave_Integrand_r = new TF1("fTsallisBlastWave_Integrand_r", TsallisBlastWave_Integrand_r, 0., 1., 8);
    if (!fTsallisBlastWave_Integrand_phi)
        fTsallisBlastWave_Integrand_phi = new TF1("fTsallisBlastWave_Integrand_phi", TsallisBlastWave_Integrand_phi, -TMath::Pi(), TMath::Pi(), 0);
    if (!fTsallisBlastWave_Integrand_y)
        fTsallisBlastWave_Integrand_y = new TF1("fTsallisBlastWave_Integrand_y", TsallisBlastWave_Integrand_y, -0.5, 0.5, 0);

    fTsallisBlastWave_Integrand_r->SetParameters(mt, pt, beta_max, temp, n, q, 0., 0.);
    Double_t integral = fTsallisBlastWave_Integrand_y->Integral(-0.5, 0.5);
    return norm * pt * integral;
}

TF1 *TsallisBlastWave(const Char_t *name, Double_t mass, Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t q = 2., Double_t norm = 1.e6)
{

    TF1 *fTsallisBlastWave = new TF1(name, TsallisBlastWave_Func, 0., 10., 6);
    fTsallisBlastWave->SetParameters(mass, beta_max, temp, n, q, norm);
    fTsallisBlastWave->SetParNames("mass", "beta_max", "T", "n", "q", "norm");
    fTsallisBlastWave->FixParameter(0, mass);
    fTsallisBlastWave->SetParLimits(1, 0.01, 0.99);
    fTsallisBlastWave->SetParLimits(2, 0.01, 1.);
    fTsallisBlastWave->SetParLimits(3, 0.1, 10.);
    fTsallisBlastWave->SetParLimits(4, 1., 10.);
    return fTsallisBlastWave;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/

TF1 *BGBlastWave_SingleFit(TH1 *h, Double_t mass, Option_t *opt = "")
{

    TF1 *f = BGBlastWave(Form("fBGBW_%s", h->GetName()), mass);
    h->Fit(f);
    h->Fit(f);
    h->Fit(f, opt);
    return f;
}

Int_t nBW;
TF1 *fBGBW[1000];
TF1 *fBGBWratio[1000];
TGraphErrors *gBW[1000];

void BGBlastWave_FCN(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

TObjArray *
BGBlastWave_GlobalFit(TObjArray *data, Double_t *mass, Double_t profile = 0.5, Bool_t computeCont = kFALSE, Bool_t fixProfile = kFALSE)
{

    /* get data */
    Int_t ndf = 0;
    nBW = data->GetEntries();
    for (Int_t idata = 0; idata < nBW; idata++)
    {
        gBW[idata] = (TGraphErrors *)data->At(idata);
        gBW[idata]->SetName(Form("gBW%d", idata));
        ndf += gBW[idata]->GetN();
    }

    /* init BG blast-wave functions */
    for (Int_t idata = 0; idata < nBW; idata++)
    {
        printf("init BG-BlastWave function #%d: mass = %f\n", idata, mass[idata]);
        fBGBW[idata] = BGBlastWave(Form("fBGBW%d", idata), mass[idata]);
    }

    if (computeCont)
        printf("-> compute contours requested\n");

    /* display data */
    TCanvas *cBW = new TCanvas("cBW");
    cBW->Divide(nBW, 1);
    for (Int_t idata = 0; idata < nBW; idata++)
    {
        cBW->cd(idata + 1);
        gBW[idata]->Draw("ap*");
    }
    cBW->Update();

    /* init minuit: nBW normalizations + 3 (beta, T, n) BG-BlastWave params */
    const Int_t nbwpars = 3;
    const Int_t nfitpars = nBW + nbwpars;
    TMinuit *minuit = new TMinuit(nfitpars);
    minuit->SetFCN(BGBlastWave_FCN);
    Double_t arglist[10];
    Int_t ierflg = 0;
    arglist[0] = 1;
    minuit->mnexcm("SET ERR", arglist, 1, ierflg);
    for (Int_t idata = 0; idata < nBW; idata++)
    {
        minuit->mnparm(idata, Form("norm%d", idata), 1.e6, 1., 0., 0., ierflg);
        ndf--;
    }
    //  minuit->mnparm(nBW + 0, "<beta>", 0.55, 0.01, 0., 1., ierflg);
    //  minuit->mnparm(nBW + 1, "T", 0.14, 0.01, 0., 1., ierflg);
    //  minuit->mnparm(nBW + 2, "n", profile, 0.1, 0., 10., ierflg);

    minuit->mnparm(nBW + 0, "<beta>", 0.7, 0.01, 0.2, 0.7, ierflg);
    minuit->mnparm(nBW + 1, "T", 0.07, 0.001, 0.07, 0.2, ierflg);
    minuit->mnparm(nBW + 2, "n", profile, 0.1, 0.6, 5., ierflg);

    ndf -= 3;
    if (fixProfile)
    {
        minuit->FixParameter(nBW + 2);
        ndf++;
    }

    /* set strategy */
    arglist[0] = 1;
    minuit->mnexcm("SET STRATEGY", arglist, 1, ierflg);

    /* start MIGRAD minimization */
    arglist[0] = 500000;
    arglist[1] = 1.;
    minuit->mnexcm("MIGRAD", arglist, 2, ierflg);

    /* set strategy */
    arglist[0] = 2;
    minuit->mnexcm("SET STRATEGY", arglist, 1, ierflg);

    /* start MIGRAD minimization */
    arglist[0] = 500000;
    arglist[1] = 1.;
    minuit->mnexcm("MIGRAD", arglist, 2, ierflg);

    /* start IMPROVE minimization */
    arglist[0] = 500000;
    minuit->mnexcm("IMPROVE", arglist, 1, ierflg);

    /* start MINOS */
    arglist[0] = 500000;
    arglist[1] = nBW + 1;
    arglist[2] = nBW + 2;
    arglist[3] = nBW + 3;
    minuit->mnexcm("MINOS", arglist, 4, ierflg);

    /* print results */
    Double_t amin, edm, errdef;
    Int_t nvpar, nparx, icstat;
    minuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
    minuit->mnprin(4, amin);

    /* get parameters */
    Double_t beta, betae, betaeplus, betaeminus, betagcc, temp, tempe, tempeplus, tempeminus, tempgcc, prof, profe, profeplus, profeminus, profgcc;
    minuit->GetParameter(nBW + 0, beta, betae);
    minuit->mnerrs(nBW + 0, betaeplus, betaeminus, betae, betagcc);
    minuit->GetParameter(nBW + 1, temp, tempe);
    minuit->mnerrs(nBW + 1, tempeplus, tempeminus, tempe, tempgcc);
    minuit->GetParameter(nBW + 2, prof, profe);
    minuit->mnerrs(nBW + 2, profeplus, profeminus, profe, profgcc);
    Double_t beta_max = 0.5 * (2. + prof) * beta;
    Double_t norm[1000], norme[1000];
    for (Int_t idata = 0; idata < nBW; idata++)
        minuit->GetParameter(idata, norm[idata], norme[idata]);

    /* printout */
    printf("[x] *********************************\n");
    printf("[x] beta_max = %f\n", beta_max);
    printf("[x] <beta>   = %f +- %f (e+ = %f, e- = %f)\n", beta, betae, betaeplus, betaeminus);
    printf("[x] T        = %f +- %f (e+ = %f, e- = %f)\n", temp, tempe, tempeplus, tempeminus);
    printf("[x] n        = %f +- %f (e+ = %f, e- = %f)\n", prof, profe, profeplus, profeminus);
    printf("[x] chi2     = %f\n", amin);
    printf("[x] ndf      = %d\n", ndf);

    /* 1-sigma contour */
    minuit->SetErrorDef(1);
    TGraph *gCont1 = NULL;
    if (computeCont)
        gCont1 = (TGraph *)minuit->Contour(50, nBW + 0, nBW + 1);
    if (gCont1)
        gCont1->SetName("gCont1");

    /* 2-sigma contour */
    minuit->SetErrorDef(4);
    TGraph *gCont2 = NULL;
    //  if (computeCont) gCont2 = (TGraph *) minuit->Contour(50, nBW + 0, nBW + 1);
    if (gCont2)
        gCont2->SetName("gCont2");

    /* display fits */
    for (Int_t idata = 0; idata < nBW; idata++)
    {
        cBW->cd(idata + 1);
        fBGBW[idata]->SetParameter(4, norm[idata]);
        fBGBW[idata]->SetParameter(1, beta_max);
        fBGBW[idata]->SetParameter(2, temp);
        fBGBW[idata]->SetParameter(3, prof);
        fBGBW[idata]->Draw("same");
    }
    cBW->Update();

    /* histo params */
    TH1D *hBW = new TH1D("hBW", "", 4, 0., 4.);
    hBW->SetBinContent(1, beta);
    hBW->SetBinError(1, betae);
    hBW->SetBinContent(2, temp);
    hBW->SetBinError(2, tempe);
    hBW->SetBinContent(3, prof);
    hBW->SetBinError(3, profe);
    hBW->SetBinContent(4, amin / ndf);

    /* BW graph */
    TGraphAsymmErrors *gBetaT = new TGraphAsymmErrors();
    gBetaT->SetName("gBetaT");
    gBetaT->SetPoint(0, beta, temp);
    gBetaT->SetPointEXlow(0, TMath::Abs(betaeminus));
    gBetaT->SetPointEXhigh(0, TMath::Abs(betaeplus));
    gBetaT->SetPointEYlow(0, TMath::Abs(tempeminus));
    gBetaT->SetPointEYhigh(0, TMath::Abs(tempeplus));

    /* prepare output array */
    TObjArray *outoa = new TObjArray();
    for (Int_t idata = 0; idata < nBW; idata++)
    {
        outoa->Add(gBW[idata]);
        outoa->Add(fBGBW[idata]);
    }
    outoa->Add(cBW);
    outoa->Add(hBW);
    outoa->Add(gBetaT);
    if (gCont1)
        outoa->Add(gCont1);
    if (gCont2)
        outoa->Add(gCont2);

    return outoa;
}

void BGBlastWave_FCN(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

    /* beta -> beta_max */
    Double_t beta = par[nBW + 0];
    Double_t T = par[nBW + 1];
    Double_t n = par[nBW + 2];
    Double_t beta_max = 0.5 * (2. + n) * beta;
#if 1
    /* check beta_max */
    if (beta_max >= 1. || beta_max <= 0.)
    {
        f = kMaxInt;
        return;
    }
    /* check T */
    if (T <= 0.)
    {
        f = kMaxInt;
        return;
    }
#endif

    Double_t pt, pte, val, vale, func, pull, chi = 0;
    /* loop over all the data */
    for (Int_t iBW = 0; iBW < nBW; iBW++)
    {
        /* set BGBW parameters */
        fBGBW[iBW]->SetParameter(4, par[iBW]);
        fBGBW[iBW]->SetParameter(1, beta_max);
        fBGBW[iBW]->SetParameter(2, T);
        fBGBW[iBW]->SetParameter(3, n);
        /* loop over all the points */
        for (Int_t ipt = 0; ipt < gBW[iBW]->GetN(); ipt++)
        {
            pt = gBW[iBW]->GetX()[ipt];
            pte = gBW[iBW]->GetEX()[ipt];
            val = gBW[iBW]->GetY()[ipt];
            vale = gBW[iBW]->GetEY()[ipt];
            func = fBGBW[iBW]->Eval(pt);
            //      func = fBGBW[iBW]->Integral(pt - pte, pt + pte);
            pull = (val - func) / vale;
            chi += pull * pull;
        }
    }

    f = chi;
}

/*****************************************************************/

#endif