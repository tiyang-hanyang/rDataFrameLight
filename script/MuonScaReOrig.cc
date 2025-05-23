#include <boost/math/special_functions/erf.hpp>
struct CrystalBall{
    double pi=3.14159;
    double sqrtPiOver2=sqrt(pi/2.0);
    double sqrt2=sqrt(2.0);
    double m;
    double s;
    double a;
    double n;
    double B;
    double C;
    double D;
    double N;
    double NA;
    double Ns;
    double NC;
    double F;
    double G;
    double k;
    double cdfMa;
    double cdfPa;
CrystalBall():m(0),s(1),a(10),n(10){
    init();
}
CrystalBall(double mean, double sigma, double alpha, double n)
    :m(mean),s(sigma),a(alpha),n(n){
    init();
}
void init(){
    double fa = fabs(a);
    double ex = exp(-fa*fa/2);
    double A  = pow(n/fa, n) * ex;
    double C1 = n/fa/(n-1) * ex; 
    double D1 = 2 * sqrtPiOver2 * erf(fa/sqrt2);
    B = n/fa-fa;
    C = (D1+2*C1)/C1;   
    D = (D1+2*C1)/2;   
    N = 1.0/s/(D1+2*C1); 
    k = 1.0/(n-1);  
    NA = N*A;       
    Ns = N*s;       
    NC = Ns*C1;     
    F = 1-fa*fa/n; 
    G = s*n/fa;    
    cdfMa = cdf(m-a*s);
    cdfPa = cdf(m+a*s);
}
double pdf(double x) const{ 
    double d=(x-m)/s;
    if(d<-a) return NA*pow(B-d, -n);
    if(d>a) return NA*pow(B+d, -n);
    return N*exp(-d*d/2);
}
double pdf(double x, double ks, double dm) const{ 
    double d=(x-m-dm)/(s*ks);
    if(d<-a) return NA/ks*pow(B-d, -n);
    if(d>a) return NA/ks*pow(B+d, -n);
    return N/ks*exp(-d*d/2);

}
double cdf(double x) const{
    double d = (x-m)/s;
    if(d<-a) return NC / pow(F-s*d/G, n-1);
    if(d>a) return NC * (C - pow(F+s*d/G, 1-n) );
    return Ns * (D - sqrtPiOver2 * erf(-d/sqrt2));
}
double invcdf(double u) const{
    if(u<cdfMa) return m + G*(F - pow(NC/u, k));
    if(u>cdfPa) return m - G*(F - pow(C-u/NC, -k) );
    return m - sqrt2 * s * boost::math::erf_inv((D - u/Ns )/sqrtPiOver2);
}
};

double get_rndm(double eta, float nL) {

    // obtain parameters from correctionlib
    double mean = cset->at("cb_params")->evaluate({abs(eta), nL, 0});
    double sigma = cset->at("cb_params")->evaluate({abs(eta), nL, 1});
    double n = cset->at("cb_params")->evaluate({abs(eta), nL, 2});
    double alpha = cset->at("cb_params")->evaluate({abs(eta), nL, 3});
    
    // instantiate CB and get random number following the CB
    CrystalBall cb(mean, sigma, alpha, n);
    TRandom3 rnd(time(0));
    double rndm = gRandom->Rndm();
    return cb.invcdf(rndm);
}


double get_std(double pt, double eta, float nL) {

    // obtain paramters from correctionlib
    double param_0 = cset->at("poly_params")->evaluate({abs(eta), nL, 0});
    double param_1 = cset->at("poly_params")->evaluate({abs(eta), nL, 1});
    double param_2 = cset->at("poly_params")->evaluate({abs(eta), nL, 2});

    // calculate value and return max(0, val)
    double sigma = param_0 + param_1 * pt + param_2 * pt*pt;
    if (sigma < 0) sigma = 0;
    return sigma; 
}


double get_k(double eta, string var) {

    // obtain parameters from correctionlib
    double k_data = cset->at("k_data")->evaluate({abs(eta), var});
    double k_mc = cset->at("k_mc")->evaluate({abs(eta), var});

    // calculate residual smearing factor
    // return 0 if smearing in MC already larger than in data
    double k = 0;
    if (k_mc < k_data) k = sqrt(k_data*k_data - k_mc*k_mc);
    return k;
}


ROOT::VecOps::RVec<float> pt_resol(ROOT::VecOps::RVec<float> pt, ROOT::VecOps::RVec<float> eta, ROOT::VecOps::RVec<unsigned char> nL)
{
    // vectorize processing
    ROOT::VecOps::RVec<float> ptResolCorr(pt.size());
    for (size_t i = 0; i < pt.size(); ++i) {
        // load correction values
        float rndm = get_rndm(eta[i], nL[i]);
        float std = get_std(pt[i], eta[i], nL[i]);
        float k = get_k(eta[i], "nom");

        // calculate corrected value and return original value if a parameter is nan
        float ptc = pt[i] * ( 1 + k * std * rndm);
        if (isnan(ptc)) ptc = pt[i];
        ptResolCorr[i] = ptc;
    }

    return ptResolCorr;
}

ROOT::VecOps::RVec<float> pt_resol_var(ROOT::VecOps::RVec<float> pt_woresol, ROOT::VecOps::RVec<float> pt_wresol, ROOT::VecOps::RVec<float> eta, string updn)
{
    // vectorize processing
    ROOT::VecOps::RVec<float> ptResolVar(pt_woresol.size());
    for (size_t i = 0; i < pt_woresol.size(); ++i) 
    {
        float k = get_k(eta[i], "nom");
        if (k == 0.0)
        {
            ptResolVar[i] = pt_wresol[i];
            continue;
        }

        float k_unc = cset->at("k_mc")->evaluate({abs(eta[i]), "stat"});
        float std_x_rndm = (pt_wresol[i] / pt_woresol[i] - 1) / k;

        float pt_var = pt_wresol[i];
        if (updn=="up"){
            pt_var = pt_woresol[i] * (1 + (k+k_unc) * std_x_rndm);
        }
        else if (updn=="dn"){
            pt_var = pt_woresol[i] * (1 + (k-k_unc) * std_x_rndm);
        }
        else {
            cout << "ERROR: updn must be 'up' or 'dn'" << endl;
        }
        ptResolVar[i] = pt_var;
    }

    return ptResolVar;
}

ROOT::VecOps::RVec<float> pt_scale(bool is_data, ROOT::VecOps::RVec<float> pt, ROOT::VecOps::RVec<float> eta, ROOT::VecOps::RVec<float> phi, ROOT::VecOps::RVec<int> charge) {
    // use right correction
    string dtmc = "mc";
    if (is_data) dtmc = "data";

    // vectorize processing
    ROOT::VecOps::RVec<float> ptScaleCorr(pt.size());
    for (size_t i = 0; i < pt.size(); ++i) {
        float a = cset->at("a_"+dtmc)->evaluate({eta[i], phi[i], "nom"});
        float m = cset->at("m_"+dtmc)->evaluate({eta[i], phi[i], "nom"});
        ptScaleCorr[i] = 1. / (m/pt[i] + charge[i] * a);
    }

    return ptScaleCorr;
}

ROOT::VecOps::RVec<float> pt_scale_var(ROOT::VecOps::RVec<float> pt, ROOT::VecOps::RVec<float> eta, ROOT::VecOps::RVec<float> phi, ROOT::VecOps::RVec<int> charge, string updn) {
    // vectorize processing
    ROOT::VecOps::RVec<float> ptScaleVar(pt.size());
    for (size_t i = 0; i < pt.size(); ++i) {
        float stat_a = cset->at("a_mc")->evaluate({eta[i], phi[i], "stat"});
        float stat_m = cset->at("m_mc")->evaluate({eta[i], phi[i], "stat"});
        float stat_rho = cset->at("m_mc")->evaluate({eta[i], phi[i], "rho_stat"});

        float unc = pt[i] * pt[i] * sqrt(stat_m * stat_m / (pt[i]*pt[i]) + stat_a*stat_a + 2*charge[i]*stat_rho*stat_m/pt[i]*stat_a);

        float pt_var = pt[i];
        if (updn=="up"){
            pt_var = pt[i] + unc;
        }
        else if (updn=="dn"){
            pt_var = pt[i] - unc;
        }
        ptScaleVar[i] = pt_var;
    }

    return ptScaleVar;
}