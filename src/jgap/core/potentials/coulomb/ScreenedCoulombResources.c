#include "ScreenedCoulombResources.h"

#if defined(__has_embed)
// Embed the data files
const unsigned char dmol_dat[] = {
#embed "../../../../../resources/dmol-screening-fit/dmol.dat"
};
const unsigned int dmol_dat_len = sizeof(dmol_dat);

const unsigned char mp2_dat[] = {
#embed "../../../../../resources/dmol-screening-fit/mp2.dat"
};
const unsigned int mp2_dat_len = sizeof(mp2_dat);
#endif
