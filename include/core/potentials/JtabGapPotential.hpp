#ifndef JTABGAPPOTENTIAL_HPP
#define JTABGAPPOTENTIAL_HPP

#include <string_view>

using namespace std;

namespace jgap {
    class JtabGapPotential : public Potential {
    public:
        JtabGapPotential();
        void save(string_view filename) const;

    };
}

#endif //JTABGAPPOTENTIAL_HPP
