#ifndef KERNEL_HPP
#define KERNEL_HPP

#include "data/BasicDataTypes.hpp"

namespace jgap {

    template<class TDescriptorData/*sparse point*/, class TIndex/*pre-calc used in multiple sparse points*/>
    class Kernel {
        public:
        virtual ~Kernel() = default;

        // TODO:
        //virtual map<TDescriptorData, > index(shared_ptr<AtomicStructure>& structure, vector<TDescriptorData> sparsePoints);
        virtual double covariance(const AtomicStructure& structure,
                                  const TIndex &indexes,
                                  const TDescriptorData &descriptor) = 0;

        virtual vector<Vector3> derivatives(const AtomicStructure& structure,
                                            const TIndex &indexes,
                                            const TDescriptorData &descriptor) = 0;

        virtual double covariance(const TDescriptorData &d1, const TDescriptorData &d2) = 0;
    };
}

#endif //KERNEL_HPP
