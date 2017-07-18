#include <iostream>
#include <vector>
using std::cout;
using std::endl;


typedef int TH1;
typedef int TString ;

class AOD_pi0
{

public:
	template  <typename T>
	void CombinatoricOfTwoInvM(std::vector<T> & collection/*reco::RecoTauPiZero */,
	                           TString typeOfCollection,
	                           TString typeOfObjects,
	                           TString typeOfConstituences,
	                           TH1 * hist_inv_m,
	                           TH1 * hist_pt = 0);

	template  <typename T>
	void CombinatoricOfTwoInvM(std::vector<T *> & collection/*reco::RecoTauPiZero */,
	                           TString typeOfCollection,
	                           TString typeOfObjects,
	                           TString typeOfConstituences,
	                           TH1 * hist_inv_m,
	                           TH1 * hist_pt = 0);
};


template  <typename T>
void AOD_pi0::CombinatoricOfTwoInvM(std::vector<T> & collection/*reco::RecoTauPiZero */,
                                    TString typeOfCollection,
                                    TString typeOfObjects,
                                    TString typeOfConstituences,
                                    TH1 * hist_inv_m,
                                    TH1 * hist_pt)
{
	cout << "Type" << endl;
}

template  <typename T>
void AOD_pi0::CombinatoricOfTwoInvM(std::vector<T *> & collection/*reco::RecoTauPiZero */,
                                    TString typeOfCollection,
                                    TString typeOfObjects,
                                    TString typeOfConstituences,
                                    TH1 * hist_inv_m,
                                    TH1 * hist_pt)
{
	cout << "Pointer" << endl;
}


int main()
{
	AOD_pi0 test;

	std::vector<int> a;
	test.CombinatoricOfTwoInvM(a, 0, 0, 0, 0);

	std::vector<int *> b;
	test.CombinatoricOfTwoInvM(b, 0, 0, 0, 0);
}

