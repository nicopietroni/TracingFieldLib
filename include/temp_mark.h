#ifndef TEMP_MARK
#define TEMP_MARK

#include <vector>

namespace Geo
{
    struct TMark
    {
        std::vector<int> ElMark;
        int GlobMark;

        void Allocate(const int &numEl)
        {
            assert(numEl >= 0);
            ElMark = std::vector<int>(numEl, GlobMark);
        }

        void AllocateIfNeeded(const int &numEl)
        {
            assert(numEl >= 0);
            if (ElMark.size()!=numEl)
                ElMark = std::vector<int>(numEl, GlobMark);
        }

        void UnMarkAll()
        {
            GlobMark++;
        }

        void Mark(const int &IndexEl)
        {
            assert(IndexEl >= 0);
            assert(IndexEl < ElMark.size());
            ElMark[IndexEl]=GlobMark;
        }

        bool IsMarked(const int &IndexEl)
        {
            assert(IndexEl >= 0);
            assert(IndexEl < ElMark.size());
            return (ElMark[IndexEl]==GlobMark);
        }

        TMark()
        {
            GlobMark = 0;
        }
    };
}

#endif
