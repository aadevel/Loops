/* *
 * The C Protein Folding Library.
 * Copyright (C) 2009 Andres Colubri.
 * Contact: andres.colubri 'AT' gmail.com
 *
 * This library was written at the Institute of Biophysical Dynamics at the University of Chicago.
 * Gordon Center for Integrated Science, W101. 929 East 57th Street, Chicago, IL 60637, USA.
 * Homepage: http://ibd.uchicago.edu/
 * 
 * Permission to use, copy, modify, and distribute this software and its
 * documentation with or without modifications and for any purpose and
 * without fee is hereby granted, provided that any copyright notices
 * appear in all copies and that both those copyright notices and this
 * permission notice appear in supporting documentation, and that the
 * names of the contributors or copyright holders not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific prior permission.
 * 
 * THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
 * CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
 * OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOFTWARE.
 * 
 */

/**
 * Discrete simulated annealing sampler based in the cooling schedule proposed by 
 * AArts and Van Laarhoven (E. Aarts & & J. Korst: "Simulated annealing and 
 * Boltzmann Machines". Wiley & Sons, 1989. Also in the chapter 4 of N. Ansari, E. Hou:
 * "Computational Intelligence for Optimization". Kluwer, 1997).
 *
 */

#ifndef __DSA_H__
#define __DSA_H__

#include "../../core/datatypes.h"
#include <stdio.h>

struct DSAParameters
{
    IntValue binsize;
    BoolValue countermoves;
    BoolValue nudge;
    BoolValue sampleomega;
    FloatValue deftemp0;
    FloatValue defenergydev0;
    IntValue updtempint;
    FloatValue updtempd;
    FloatValue logupdtempd;
    FloatValue energyrmsdsmooth;
    BoolValue checktempconv;
    FloatValue convepsilon;
    IndexValue maxnumsteps;
    BoolValue randominit;
};

ErrorCode reset_status_parameters(void);


ErrorCode set_system_for_dsa(struct AtomData *_atoms, FloatValue *_coords, 
                                                      FloatValue *_distances, FloatValue *_invdistances, 
                                                      BoolValue *_atommask, BoolValue *_distmask, IndexValue _natoms, 
                             struct ResidueData *_residues, IndexValue _nres, 
                             struct ChainData *_chains, IndexValue _nchains);

ErrorCode init_dsa(const char* cfg);
ErrorCode finish_dsa(void);

ErrorCode pre_dsa_step(void);
ErrorCode post_dsa_step(void);

ErrorCode run_dsa_step(void);

ErrorCode load_phi_psi_dist_from_file(const char* file,FloatValue * array);
ErrorCode load_rama_map_from_file(const char* file,IndexValue *indexarray,FloatValue *maparray,BoolValue loadindex);
ErrorCode load_omega_probs_from_file(const char* file,FloatValue *omega);

ErrorCode write_current_conf(FILE* path);

#endif
