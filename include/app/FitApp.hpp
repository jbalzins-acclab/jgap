#ifndef JGAP_FITAPP_HPP
#define JGAP_FITAPP_HPP

namespace jgap {

    shared_ptr<Fit> fit(string );

    void printErrors(const string& potFileName, const string& testFile, const string& groupedBy,
                     map<string, vector<double>>& energyDifferences,
                     map<string, vector<double>>& forceDifferences,
                     map<string, vector<array<jgap::Vector3, 3>>>& virialDifferences);
    double rms(const vector<double>&);
    double deviation(const vector<double>&);
}

#endif
