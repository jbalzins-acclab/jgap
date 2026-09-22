#ifndef JGAP_SCREENEDCOULOMBRESOURCES_H
#define JGAP_SCREENEDCOULOMBRESOURCES_H

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__has_embed)
// External declarations for the embedded data
extern const unsigned char dmol_dat[];
extern const unsigned int dmol_dat_len;

extern const unsigned char mp2_dat[];
extern const unsigned int mp2_dat_len;
#else
// Provide empty fallbacks if #embed is not supported
static const unsigned char dmol_dat[] = {0};
static const unsigned int dmol_dat_len = 0;

static const unsigned char mp2_dat[] = {0};
static const unsigned int mp2_dat_len = 0;
#endif

#ifdef __cplusplus
}
#endif

#endif
