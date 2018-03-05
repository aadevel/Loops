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
 * Authors: Aashish Adhikari (Multiple chain implementation)
 *
 */



#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../../core/memory.h"
#include "../../core/modules.h"
#include "../../core/fileutils.h"
#include "../../core/errorhandler.h"
#include "../../core/rotation.h"
#include "../../core/distance.h"
#include "../../core/random.h"
#include "../../core/strutils.h"
#include "dsa.h"

#include "../shared/dsadope.h"   // Include declaring variables and constants shared between DSA and DOPE modules.
extern BoolValue *flags;         // Flags array shared with the DOPE module (has to be declared extern to be shared).
extern FloatValue T;             // Temperature variable shared with the DOPE module (has to be declared extern to be shared).

#ifndef min
#define min(x, y)	((x) < (y) ? (x) : (y))
#endif
#ifndef max
#define	max(x, y)	((x) > (y) ? (x) : (y))
#endif

// Local copy of system variables.
struct AtomData *atoms;          // Array of atom data.
FloatValue *coords;              // Array of atom coordinates.
FloatValue *distances;           // Array of atom pair-wise distances.
FloatValue *invdistances;        // Array of atom pair-wise 1/distances.
BoolValue *atommask;             // Array of atom masking values.
BoolValue *distmask;             // Array of distance masking values.  
IndexValue natoms;               // Total number of atoms.

struct ResidueData *residues;    // Array of residue data.
IndexValue nres;                 // Total number of residues.
 
struct ChainData *chains;        // Array of chain data.
IndexValue nchains;              // Total number of chains.

static const IndexValue maxnumenergymods=50;

time_t starttime;

IndexValue *energymods;
FloatValue *energycoeffs;
IndexValue nenergymods;
FloatValue nudgenoise;

//IndexValue rbasinsmod;
IndexValue roundnum;

FloatValue phi_psi[2];
FloatValue psi_phi_store[2];
FloatValue *rama_distribution;
IndexValue *rama_index;
FloatValue *omega_probs;

FloatValue *psi_phi;
FloatValue *psi_psi;
FloatValue *phi_phi;
FloatValue *phi_psi_long;

IndexValue totalnumbins;
IndexValue totalpsiphibins;


//aashish
//BoolValue multchain;
//IndexValue 
BoolValue ligatechains;

//BoolValue rangeset,restosample;
//IndexValue rangemin=144;
//IndexValue rangemax=154;

//ignore this
//int restomovearray[]={22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49 };

//int restomovearray[]={7,8,9,10,11,12,13,14,15,16,17,18,19, 39,40,41,42,43,44,45,46,47,48 };

//int restomovearray[]={44,45,46,54,55,56,57,58,59,60,61,62,63,64,65};

//int restomovearray[]={81,82,83};


//int restomovearray[]={53,54,55,56,57,58,59,60};
int *restomovearray = NULL;
int nrestosample = 0;
BoolValue restosample;

//IndexValue nrestomove1=0;
//IndexValue nrestomove2=6;


BoolValue writetrajectory, saveminimal, saveoptimal;
char* trajectoryname;
BoolValue convergence,finished;
IntValue dsastep;
IntValue ndecsteps, nincsteps, nrejsteps;
IntValue nramachanges;
FloatValue temp, prevtemp; 
FloatValue curenergy;
FloatValue *curenergies;
FloatValue optenergy;
FloatValue energysdevsmooth;
IntValue nsolutions;
FloatValue energysum;
FloatValue energysumsq;

FloatValue *curconf;
FloatValue *optconf;

#define rama_index_pos(res, bin) ((res)* totalnumbins + bin)
#define psi_phi_index_pos(aa0, aa1, bin) (aa0 + 20 * aa1 + 400 * bin)
#define aa_index_pos(aa0, aa1) (aa0 + 20 * aa1)


typedef struct cNode{
  char * m_Token;
  struct cNode * m_pNext;
} CNode;


struct DSAParameters followuppar;

static struct DSAParameters par =
{
    5,       // binsize
    FALSE,   // countermoves
    FALSE,    // nudge
    FALSE,   // sampleomega
    100.0,   // deftemp0
    100,       // defenergydev0
    100,     // updtempint
    0.1,     // updtempd
    0.0,     // logupdtempd
    0.8,     // energyrmsdsmooth;
    TRUE,    // checktempconv;
    1,    // convepsilon
    0,        //maxnumsteps
    TRUE        //randominit
};

ErrorCode save_current_conf(IndexValue res, IndexValue chn)
{
    IndexValue firstat;
    
    if(par.countermoves || par.nudge){firstat = residues[res-1].firstatom;}
    else{ firstat = residues[res].bonds[PHI].atom2;}

    // The last atom affected by the rotation should be always the last atom in the chain.
    IndexValue lastat = residues[chains[chn].lastres].lastatom;

    IndexValue nat = 1 + lastat - firstat;
    //printf("res %i firstat %i lastat %i nat %i\n",res,firstat,lastat,nat);
    memcpy(optconf + 4 * firstat, coords + 4 * firstat, 4 * sizeof(FloatValue) * nat);
    return NO_ERROR;
}

ErrorCode restore_prev_conf(IndexValue res, IndexValue chn)
{
    IndexValue firstat;

    if(par.countermoves || par.nudge){firstat = residues[res-1].firstatom;}
    else{ firstat = residues[res].bonds[PHI].atom2;}

    // The last atom affected by the rotation should be always the last atom in the chain.
    IndexValue lastat = residues[chains[chn].lastres].lastatom;
    IndexValue nat = 1 + lastat - firstat;
    //printf("res %i firstat %i lastat %i nat %i\n",res,firstat,lastat,nat);
    memcpy(coords + 4 * firstat, optconf + 4 * firstat, 4 * sizeof(FloatValue) * nat);    
    
     //aashish
      return calculate_masked_distances(distances, coords, natoms, distmask);
     // return calculate_all_distances(distances,coords,natoms);
}

ErrorCode save_opt_conf()
{
    memcpy(optconf, coords, 4 * sizeof(FloatValue) * natoms);
    return NO_ERROR;
}

ErrorCode restore_opt_conf()
{
    memcpy(coords, optconf, 4 * sizeof(FloatValue) * natoms);
    return calculate_all_distances(distances, coords, natoms);
}

void choose_psi_phi(IndexValue res)
{
    FloatValue rama_prob, ref_probs;
    IndexValue i;

    rama_prob = get_frandom();
    ref_probs = 0.0;
    IndexValue aa = residues[res].type - 1;
    IndexValue aa1 = residues[res + 1].type - 1;
    for (i = 0; i < 324; i++)
    {
	    if ((rama_prob >= ref_probs) && (rama_prob < psi_phi[psi_phi_index_pos(aa, aa1, i)] + ref_probs))
	    {
	        psi_phi_store[0] = (20.0 * (i % 18)) - 180.00 + (20.00 * get_frandom()); 
	        psi_phi_store[1] = (20.0 * ((i - i % 18) / 18)) - 180.00 + (20.00 * get_frandom());
	        break;
	    }
	    else ref_probs += psi_phi[psi_phi_index_pos(aa, aa1, i)];
    }
    //printf("Choosing phi, psi: res: %i phi: %f psi: %f rama_prob: %f\n",res,phi_psi[0],phi_psi[1],rama_prob);
}

void choose_phi_psi(IndexValue res)
{
    FloatValue rama_prob, ref_probs;
    IndexValue i, bin;

    rama_prob = get_frandom();
    ref_probs = 0.0;
    for (i = 0; i < totalnumbins; i++)
    {
	    if ((rama_prob >= ref_probs) && (rama_prob < rama_distribution[rama_index_pos(res, i)] + ref_probs))
	    {
            bin = rama_index[rama_index_pos(res, i)];

	        phi_psi[0] = (5.0 * (bin % 72)) - 180.00 + (5.00 * get_frandom()); 
	        phi_psi[1] = (5.0 * ((bin - bin % 72) / 72)) - 180.00 + (5.00 * get_frandom());
	        break;
        }
        else ref_probs += rama_distribution[rama_index_pos(res, i)];
    }
    //printf("Choosing phi, psi: res: %i phi: %f psi: %f rama_prob: %f\n",res,phi_psi[0],phi_psi[1],rama_prob);
}

IndexValue general_bin(FloatValue ang1, FloatValue ang2)
{
    IndexValue ang1_bin, ang2_bin, bin;
    ang1_bin = (IndexValue)(ang1 + 180.0) / 20.0; // Max value: 17.
    ang2_bin = (IndexValue)(ang2 + 180.0) / 20.0; // Max value: 17.
    bin = ang1_bin + 18 * ang2_bin; // Max value: 18 * 18 - 1
    return bin;
}

void set_phi_psi(IndexValue res, FloatValue phi, FloatValue psi)
{
    residues[res].bonds[PHI].newangle = phi;
    residues[res].bonds[PSI].newangle = psi;
}

void set_psi_phi(IndexValue res, FloatValue psi, FloatValue phi)
{
    residues[res].bonds[PSI].newangle = psi;
    residues[res + 1].bonds[PHI].newangle = phi;
}

ErrorCode generate_rand_conf()
{
    IndexValue res, count;
    FloatValue phi_left, psi_left, old_phi, old_psi, burial, burial_new;
    BoolValue psi_phi_ok, phi_phi_ok, psi_psi_ok, phi_psi_long_ok, second_pos, last_pos, donehere;
    int chn;

    donehere = FALSE;

    for ( chn=0; chn < nchains; chn++){
   
   	while (TRUE)
    	{
	
//	printf("Chose chn %i \n" , chn);
	
	choose_psi_phi(chains[chn].firstres);
        set_psi_phi(chains[chn].firstres, psi_phi_store[0], psi_phi_store[1]);

  //      printf("%f %f \n",psi_phi_store[0], psi_phi_store[1]);
        
        for (res = chains[chn].firstres + 1; res <= chains[chn].lastres; res++)
        {
	        second_pos = (res == chains[chn].firstres + 1);
	        last_pos = (res == chains[chn].lastres);
	        count = 0;
	        
	        while (TRUE)
	        {
		choose_phi_psi(res);
                IndexValue aa0 = residues[res - 1].type - 1;
                IndexValue aa = residues[res].type - 1;

	            psi_left = residues[res - 1].bonds[PSI].angle;
		        phi_left = residues[res - 1].bonds[PHI].angle;
	            psi_phi_ok = (psi_phi[psi_phi_index_pos(aa0, aa, general_bin(psi_left, phi_psi[0]))] > 0);
		        psi_psi_ok = (psi_psi[psi_phi_index_pos(aa0, aa, general_bin(psi_left, phi_psi[1]))] > 0);
		        phi_phi_ok = (phi_phi[psi_phi_index_pos(aa0, aa, general_bin(phi_left, phi_psi[0]))] > 0);
		        phi_psi_long_ok = (phi_psi_long[psi_phi_index_pos(aa0, aa, general_bin(phi_left, phi_psi[1]))] > 0);
		        if ((second_pos && psi_phi_ok && psi_psi_ok) || 
                    (last_pos && psi_phi_ok && phi_phi_ok) || 
                    (!last_pos && !second_pos && psi_phi_ok && psi_psi_ok && phi_phi_ok && phi_psi_long_ok))
	            {
		            set_phi_psi(res, phi_psi[0], phi_psi[1]);
		            donehere = last_pos;
	                break;
	            }
		        else count++;
                
		        if (count == 1000) break;
	        }
            
	        if (count == 1000) break;
        }
	    if (donehere) break;
    }
    }

    update_all_coords(coords, distmask, atoms, natoms, residues, nres, chains, nchains, TRUE);
    calculate_all_distances(distances, coords, natoms);
    return NO_ERROR;
}


BoolValue change_rama(IndexValue res, IndexValue chn)
{
    FloatValue new_phi, new_psi;
    FloatValue phi_L, psi_L, phi_R, psi_R;
    FloatValue old_phi, old_psi, delta_phi, delta_psi;
    FloatValue noise_phi, noise_psi, noise_omega;
    FloatValue omega, omega_prob, omega_1;
    FloatValue pep_psi, pep_phi, left_phi, left_psi, right_phi, right_psi;
    BoolValue omega_cis, final_pos, first_pos, second_pos, penult_pos;
    BoolValue L_psi_phi_ok, L_phi_phi_ok, L_psi_psi_ok, L_phi_psi_long_ok; 
    BoolValue R_psi_phi_ok, R_phi_phi_ok, R_psi_psi_ok, R_phi_psi_long_ok;

    omega_cis = FALSE;
    nramachanges++;

    //printf("Init Change rama");

    save_current_conf(res, chn); 
    clear_distmask(distmask, natoms);

   // printf(" chains is %i and res is %i \n", chn, res);

    
    if (!par.countermoves && !par.nudge)
    {
        first_pos = (res == chains[chn].firstres);
        second_pos = (res == chains[chn].firstres + 1);
        final_pos = (res == chains[chn].lastres);
        penult_pos = (res == chains[chn].lastres - 1);

        if (!first_pos)
        {
            left_psi = residues[res - 1].bonds[PSI].angle;
            left_phi = residues[res - 1].bonds[PHI].angle;
        }
        if (!final_pos)
        {
            right_psi = residues[res + 1].bonds[PSI].angle;
            right_phi = residues[res + 1].bonds[PHI].angle;
        }

        IntValue move_count = 0;
        while (TRUE)
        {
            choose_phi_psi(res);
            move_count += 1;
 
            if (move_count == 1000) break; // Ok, we didn't find any phi, psi pair that gives an allowed psi, phi pair, so we'll just take it anyway

            IndexValue aa = residues[res].type - 1;           
            if (!first_pos)
            {
                IndexValue aa0 = residues[res - 1].type - 1;
                L_psi_phi_ok = (psi_phi[psi_phi_index_pos(aa0, aa, general_bin(left_psi, phi_psi[0]))] > 0);
	        L_psi_psi_ok = (psi_psi[psi_phi_index_pos(aa0, aa, general_bin(left_psi, phi_psi[1]))] > 0);
	        L_phi_phi_ok = (phi_phi[psi_phi_index_pos(aa0, aa, general_bin(left_phi, phi_psi[0]))] > 0);
                L_phi_psi_long_ok = (phi_psi_long[psi_phi_index_pos(aa0, aa, general_bin(left_phi, phi_psi[1]))] > 0);
	    }

            if (!final_pos)
            {
                IndexValue aa1 = residues[res + 1].type - 1;
	        R_psi_phi_ok = (psi_phi[psi_phi_index_pos(aa, aa1, general_bin(phi_psi[1], right_phi))] > 0);
	        R_psi_psi_ok = (psi_psi[psi_phi_index_pos(aa, aa1, general_bin(phi_psi[1], right_psi))] > 0);
	        R_phi_phi_ok = (phi_phi[psi_phi_index_pos(aa, aa1, general_bin(phi_psi[0], right_phi))] > 0);
	        R_phi_psi_long_ok = (phi_psi_long[psi_phi_index_pos(aa, aa1, general_bin(phi_psi[0], right_psi))] > 0);
	    }

	    if (first_pos)
	    {
	        if (R_psi_phi_ok && R_psi_psi_ok) break;
	        else continue;
	    }
	    if (second_pos)
	    {
                if (L_psi_phi_ok && L_psi_psi_ok && R_psi_phi_ok && R_psi_psi_ok && R_phi_phi_ok && R_phi_psi_long_ok) break;
	        else continue;
	    }
	    if (penult_pos)
	    {
	        if (R_psi_phi_ok && R_phi_phi_ok && L_psi_phi_ok && L_psi_psi_ok && L_phi_phi_ok && L_phi_psi_long_ok) break;
	        else continue;
	    }
	    if (final_pos)
	    {
	        if (L_psi_phi_ok && L_phi_phi_ok) break;
	        else continue;
	    }
	    if (R_psi_phi_ok && R_psi_psi_ok && R_phi_phi_ok && R_phi_psi_long_ok && L_psi_phi_ok && L_psi_psi_ok && L_phi_phi_ok && L_phi_psi_long_ok) break;
	    else continue;
        }

	set_phi_psi(res, phi_psi[0], phi_psi[1]);
       //printf("Setting phi psi for res %i as %f, %f \n", res,phi_psi[0], phi_psi[1]);

	    if (par.sampleomega && !first_pos)
	    {
            IndexValue aa0 = residues[res - 1].type - 1;
            IndexValue aa = residues[res].type - 1;

	        omega_prob = get_frandom();
	        if (omega_prob < omega_probs[aa_index_pos(aa0, aa)]) omega = 0.0;
	        else omega = 180.0;

            residues[res - 1].bonds[OMEGA].newangle = omega;
	    }
    }
    else if(par.countermoves && !par.nudge)
    {
	choose_phi_psi(res);
	new_phi = phi_psi[0];
	new_psi = phi_psi[1];
   
	old_phi = residues[res].bonds[PHI].angle;
	old_psi = residues[res].bonds[PSI].angle;
	delta_phi = new_phi - old_phi;
        delta_psi = new_psi - old_psi;

        // Left compensation
        phi_L = residues[res - 1].bonds[PHI].angle;
        psi_L = residues[res - 1].bonds[PSI].angle;
        psi_L -= delta_phi;
        if (psi_L < -180) psi_L = 360 + psi_L;
        if (psi_L > 180) psi_L = psi_L - 360;
        set_phi_psi(res - 1, phi_L, psi_L);

        // Central insertion
        set_phi_psi(res, new_phi, new_psi);

        // Right compensation
        phi_R = residues[res + 1].bonds[PHI].angle;
        psi_R = residues[res + 1].bonds[PSI].angle;
        phi_R -= delta_psi;
        if (phi_R < -180) phi_R = 360 + phi_R;
        if (phi_R > 180) phi_R = phi_R - 360;
        set_phi_psi(res + 1, phi_R, psi_R);
    }
    else{
       
       old_phi = residues[res].bonds[PHI].angle;
       old_psi = residues[res].bonds[PSI].angle;

        //FloatValue noise = 1.0;
 	new_phi = old_phi + (get_frandom() - 0.5) * nudgenoise;
	new_psi = old_psi + (get_frandom() - 0.5) * nudgenoise;
	set_phi_psi(res, new_phi, new_psi);


    }


    update_all_coords(coords, distmask, atoms, natoms, residues, nres, chains, nchains, TRUE);
    calculate_masked_distances(distances, coords, natoms, distmask);

    update_all_torsional_angles(residues, chains, nchains, coords);

    return TRUE;
}

ErrorCode generate_rand_selected()
{
   IndexValue resi, chni;
   
   int pick;

   for ( pick=0; pick < nrestosample; pick++){
 
//   printf("Pick is %i \n", pick);

  resi = restomovearray[pick] - 1;    
  chni = residues[resi].chain;


 //  printf("Changing rama for residue %i \n", resi);
   change_rama(resi, chni);

   }

    update_all_coords(coords, distmask, atoms, natoms, residues, nres, chains, nchains, TRUE);
    calculate_all_distances(distances, coords, natoms);


   return NO_ERROR;
}


BoolValue generate_new_conf(IndexValue *res, IndexValue *chn)
{
    IndexValue res0, res1;

    // This code only works for the first chain (0).
//    res0 = chains[*chn].firstres;
 //   res1 = chains[*chn].lastres; 

      // if(rangeset){
//	res0= rangemin -1;
//	res1= rangemax -1;
  //     }
  //     else{
       res0 = 0;
       res1 = nres - 1;
//	}
    if (par.countermoves || par.nudge)
    {
        res0++;
        res1--;
    }

   if(restosample){
   int picked = get_irandom(0,nrestosample-1);
   *res = restomovearray[picked] -1 ;
   }
   else{
    *res = get_irandom(res0 + 1, res1 - 1);
//    *chn = get_irandom(0,nchains-1);
    }
    *chn = residues[*res].chain;

    //printf(" zone3 %i %i \n", *res, *chn);

    return change_rama(*res, *chn);
}

ErrorCode calculate_mean_and_sdev(FloatValue *energymean, FloatValue *energysdev)
{
    if (0 < nsolutions)
    {
        *energymean = energysum / nsolutions;
        if (1 < nsolutions)
        {
            *energysdev = (energysumsq - nsolutions * (*energymean) * (*energymean)) / (nsolutions - 1);
            *energysdev = 0 < *energysdev ? sqrt(*energysdev) : 0;
        }
        else *energysdev = 0;
    }
    else
    {
        *energymean = 0;
        *energysdev = 0;
    }

   return NO_ERROR;
}

ErrorCode set_system_for_dsa(struct AtomData *_atoms, FloatValue *_coords, 
                                                      FloatValue *_distances, FloatValue *_invdistances, 
                                                      BoolValue *_atommask, BoolValue *_distmask, IndexValue _natoms, 
                             struct ResidueData *_residues, IndexValue _nres, 
                             struct ChainData *_chains, IndexValue _nchains)
{
    atoms = _atoms;
    coords = _coords;
    distances = _distances;
    invdistances = _invdistances;
    atommask = _atommask;
    distmask = _distmask;
    natoms = _natoms;
    residues = _residues;
    nres = _nres;
    chains = _chains;
    nchains = _nchains;

    return NO_ERROR;
}

ErrorCode load_omega_probs_from_file(const char* file,FloatValue *omega){
    if (load_text_file_into_buffer(file) == NO_ERROR)
    {
        printf("        Loading %s...",file);
        while (read_input_line_from_buffer())
        {
             IndexValue aa0,aa1;
             FloatValue prob;
             get_input_ival(&aa0," ", 0); 
             get_input_ival(&aa1," ", 1); 
             get_input_fval(&prob," ",2);
             omega[aa_index_pos(aa0,aa1)]=prob;
        }
        delete_input_buffer();
        printf("ok.\n");
    }
    return NO_ERROR;
}

ErrorCode load_rama_map_from_file(const char* file,IndexValue *indexarray,FloatValue *maparray,BoolValue loadindex){
    IndexValue currentres,bin;
    IndexValue currentivalue,count,currentindex;
    FloatValue currentfvalue;
    IndexValue initmaxlinelength=get_max_line_length();

    set_max_line_length(1000000);
    //printf("using totalnumbins: %i\n",totalnumbins);
    if (load_text_file_into_buffer(file) == NO_ERROR)
    {
        printf("        Loading %s...",file);
        currentres=0;
        //printf("        loadindex: %i  indexarray!=NULL %i maparray!=NULL: %i\n",loadindex,indexarray!=NULL,maparray!=NULL);
        while (read_input_line_from_buffer())
        {
            for(count=0;count<totalnumbins;count++){
                //bin=(count%totalnumbins);
                //currentres=count/totalnumbins;
                if(loadindex && (indexarray != NULL))
                {
                    get_input_ival(&currentivalue," ", count);
                    currentindex=rama_index_pos(currentres,count);
                    indexarray[currentindex]=currentivalue;
                }
                else if (!loadindex && (maparray !=NULL) ) 
                {
                    get_input_fval(&currentfvalue, " ", count);
                    currentindex=rama_index_pos(currentres,count);
                    maparray[currentindex]=currentfvalue;
                }
                else 
                {
                    print_error(ITEM_NOT_FOUND_ERROR,"load_rama_map_from_file",0);
                    return ITEM_NOT_FOUND_ERROR;
                }
            }
            currentres++;
        }
        delete_input_buffer();
        printf("ok.\n");
    }
    set_max_line_length(initmaxlinelength);
    return NO_ERROR;
}


ErrorCode load_phi_psi_dist_from_file(const char* file,FloatValue * array){
    if (load_text_file_into_buffer(file) == NO_ERROR)
    {
        printf("        Loading %s...",file);
        IndexValue tocount=18*18+2,index,currentindex=0;
        IndexValue indices[2];
        IndexValue count=0;
        while (read_input_line_from_buffer())
        {
            FloatValue currentvalue;
            IndexValue modvalue=count % tocount;
            if (modvalue==0||modvalue==1){
                get_input_ival(&index, " ", 0);
                indices[modvalue]=index;
            }
            else{
                get_input_fval(&currentvalue, " ", 0);
                //from above: #define psi_phi_index_pos(aa0, aa1, bin) (aa0 + 20 * aa1 + 400 * bin)
                currentindex=psi_phi_index_pos(indices[0],indices[1],(modvalue-2) ); // WARNING: parenthases here very important because a macro is used, not a function
                array[currentindex]=currentvalue;
//                printf("%i %i %i %i %e\n",indices[0],indices[1],modvalue-2,currentindex,currentvalue);
            }
            count++;
        }
        delete_input_buffer();
        printf("ok.\n");
    }
    return NO_ERROR;
}

ErrorCode reset_status_parameters(void){
    dsastep = 0;
    nsolutions = 0;
    energysum = 0;
    energysumsq = 0;
    convergence = FALSE;
    finished = FALSE;

    ndecsteps = 0;
    nincsteps = 0;
    nrejsteps = 0;
    nramachanges = 0;
    return NO_ERROR;
}

ErrorCode read_range(){



return NO_ERROR;
}


ErrorCode init_dsa(const char* cfg)
{
    starttime=time(NULL);
    roundnum=0;

    restosample=FALSE;
    
    nudgenoise = 5.0;

    //setting default followup par settings
    followuppar=par;

    IndexValue currentindex;
    FloatValue energycoeff;

    ErrorCode error = NO_ERROR;

    char **energycfgs = NULL;

    char *linehead = NULL;
    char *energyname = NULL;
    //char *energycfg = NULL;
//    char *rbasinsname = NULL;
//    char *rbasinscfg = NULL;

    char *ramadir = NULL;
    char *name = NULL;
    char *pardir = NULL;
    char *currentstr = NULL;

    char *seedfile=NULL;
    BoolValue loadseed=FALSE;

    trajectoryname = NULL;
    writetrajectory = 0;
    saveoptimal=1;
    saveminimal=0;

    int counter=0;


    create_index_array(&energymods,maxnumenergymods);
    create_float_array(&energycoeffs,maxnumenergymods);
    create_float_array(&curenergies,maxnumenergymods);
    nenergymods=0;

    
    
    IndexValue maxnumstepsin,i;


    printf("    Initializing DSA sampler...\n");

    if (open_input_text_file(cfg) == NO_ERROR)
    {
        while (read_input_line_from_file())
        {
            get_input_strval(&linehead, "=", 0);
            char firstchar=linehead[0];
            if (firstchar=='#') {
                free(linehead);
                linehead = NULL;
                continue;
            }
            if (strcmp(linehead, "ENERGY FUNCTION") == 0)
            {
                if (nenergymods<maxnumenergymods){
                    get_input_strval(&energyname, "=,", 1);
                    error = get_energy_module(&energymods[nenergymods], energyname);
                    if (error != NO_ERROR)
                    {
                        print_error(error, "init_dsa", 1, " when getting energy module");
                    }

                    char *energycfg = NULL;
                    // The energycfg pointer is allocated in get_input_strval and subsequently stored in energycfgs[nenergymods]. So it only needs to be de-allocated
                    // when freeing the array energycfgs.
                    get_input_strval(&energycfg,"=,",2);
                    energycfgs = (char **) realloc (energycfgs,(nenergymods+1)*sizeof(char *));                    
                    //energycfgs[nenergymods]=strdup(energycfg);
                    energycfgs[nenergymods]=energycfg;

                    IndexValue result = get_input_fval(&energycoeff,"=,",3);
                    if (result==ITEM_NOT_FOUND_ERROR) energycoeff=1.0;
                    energycoeffs[nenergymods]=energycoeff; 
                    printf("        Loading Energy Function %s with coefficient %f...\n",energyname,energycoeff);
                    nenergymods++;
                }
                else {
                    printf("        Error: No more energy mods can be added, exceeded max (%i)\n",maxnumenergymods);
                }
            }
	    
	    //aashish
	    else if (strcmp(linehead, "SAMPLE RANGE") == 0)
	    {
		restosample = TRUE;
//		get_input_strval();  
	      char * sline = NULL;
              get_input_strval(&sline, "=", 1);

	      CNode * pList = NULL;
	      char * sToken = strtok(sline, ",");

	      while( NULL != sToken){
		      
		      CNode * pNode = (CNode *) malloc(sizeof(CNode));
		      pNode->m_Token = strdup(sToken);
		      pNode->m_pNext = pList;
		      pList = pNode;

		      sToken = strtok(NULL,",");
	      }

              CNode * pNode = pList;

	      while(NULL != pNode){
              
		      if(strstr(pNode->m_Token, "-") != NULL){
		          //printf("X \n");

			  char * ssToken = strtok(pNode->m_Token, "-");
 			  int ss[2],i=0;
			  while (NULL != ssToken){
				ss[i]=atoi(strdup(ssToken));
				i++;
				ssToken = strtok(NULL,"-");
			  }
                      
		      int j=0;
		      for(j=ss[0]; j <=ss[1]; j++){
			 restomovearray = (int *) realloc (restomovearray, (counter+1)*sizeof(int));
			 restomovearray[counter++] = j;
		      }

		      }
		      else{
		  	int id=atoi(pNode->m_Token);
			// printf(" i is %i \n", id);
		       restomovearray = (int *) realloc (restomovearray, (counter+1)*sizeof(int));
		       restomovearray[counter++] = id; 
		      }


		nrestosample=counter;
	 	printf("Nrestosample is %i \n",nrestosample);
		
	      free(pNode->m_Token);
              CNode * pNext = pNode->m_pNext;
	      free(pNode);
	      pNode = pNext;
	      }

	    }
	    else if (strcmp(linehead, "NUDGE") == 0){
		     get_input_ival(&par.nudge, "=", 1);
	    }
	    else if (strcmp(linehead, "NUDGE NOISE") ==0 ){
		     get_input_ival(&nudgenoise, "=", 1);
	    }
/*
            else if (strcmp(linehead,"ENERGY CFG")==0)
            {
                get_input_strval(&energycfg, "=", 1);
            }
            else if (strcmp(linehead, "RB CALCULATOR") == 0)
            {
                get_input_strval(&rbasinsname, "=", 1);
            }
            else if (strcmp(linehead, "RB CONFIG") == 0)
            {
                get_input_strval(&rbasinscfg, "=", 1);
            }
*/
            else if (strcmp(linehead, "MAXIMUM NUMBER OF STEPS") == 0)
            {
                get_input_ival(&maxnumstepsin, "=", 1);
                par.maxnumsteps=maxnumstepsin;
            }
            else if (strcmp(linehead, "TRAJECTORY FILE") == 0)
            {
	        writetrajectory = 1;
                get_input_strval(&trajectoryname, "=", 1);
            //aashish - rand traj
	    
	      char buff [50];
	      srand(time(0));
	      int tempo = rand();
	      printf("tempo is %i" , tempo);
	      sprintf(buff, "-%i.pdt",tempo);
	     append_string(&trajectoryname, name, buff);
	    }
	    else if (strcmp(linehead,"SAVE OPTIMAL") == 0){
		saveoptimal = 1; 
	    }
	    else if (strcmp(linehead, "SAVE MINIMAL") == 0){
		saveminimal = 1;
	    }
            else if (strcmp(linehead, "PARAMETER DIR") == 0)
            {
                get_input_strval(&pardir, "=", 1);
            }
            else if (strcmp(linehead, "RANDOM INIT") == 0)
            {
                get_input_ival(&currentindex, "=", 1);
                par.randominit=currentindex;
            }
            else if (strcmp(linehead, "DEFAULT INITIAL TEMPERATURE") == 0)
            {
                get_input_ival(&currentindex, "=", 1);
                par.deftemp0=currentindex;
            }
            else if (strcmp(linehead, "TEMPERATURE UPDATE INTERVAL") == 0)
            {
                get_input_ival(&currentindex, "=", 1);
                par.updtempint=currentindex;
            }
            else if (strcmp(linehead, "RAMA DIR") == 0)
            {
                get_input_strval(&ramadir, "=", 1);
            }
            else if (strcmp(linehead, "NAME") == 0)
            {
                get_input_strval(&name, "=", 1);
            }
            else if (strcmp(linehead, "BIN SIZE") == 0)
            {
                get_input_ival(&par.binsize, "=", 1);
            }
            else if (strcmp(linehead, "COUNTER MOVES") == 0)
            {
                get_input_ival(&par.countermoves, "=", 1);
            }
            else if (strcmp(linehead, "SAMPLE OMEGA") == 0)
            {
                get_input_ival(&par.sampleomega, "=", 1);
            }
            else if (strcmp(linehead, "DEFAULT ENERGY DEVIATION") == 0)
            {
                get_input_fval(&par.defenergydev0, "=", 1);
            }
            else if (strcmp(linehead, "UPDATE TEMP COEFFICIENT") == 0)
            {
                get_input_fval(&par.updtempd, "=", 1); //logupdtempd set below
            }
            else if (strcmp(linehead, "CHECK TEMP CONVERGENCE") == 0)
            {
                get_input_ival(&par.checktempconv, "=", 1);
            }
            else if (strcmp(linehead, "SMOOTH DEVIATION COEFFICIENT") == 0)
            {
                get_input_fval(&par.energyrmsdsmooth, "=", 1);
            }
            else if (strcmp(linehead, "TEMP CONVERGENCE COEFFICIENT") == 0)
            {
                get_input_fval(&par.convepsilon, "=", 1);
            }
            else if (strcmp(linehead, "SEED FILE") == 0)
            {
                get_input_strval(&seedfile,"=",1);
                loadseed=TRUE;
            }
            else
            {
                print_error(UNKNOWN_PARAMETER_ERROR, "init_oops", 2, "header ", linehead);
            }
            free(linehead); 
            linehead = NULL; // Make sure of set this string pointer to NULL, beacause it will be used again when reading the next line. 
        }
        close_input_file();    
    }

    int xx;
    for(xx=0; xx < counter; xx++){
		printf("The element is %i \n", restomovearray[xx]);
	    }





    printf("check this");

    for(i=0;i<nenergymods;i++) execute_init_func_energy_module(energymods[i], energycfgs[i]);

    free(energyname);
    //free(energycfg);

    printf("what what");

    if (pardir==NULL) copy_string(&pardir,".");
    if (ramadir==NULL) copy_string(&ramadir,".");
    if (name==NULL) copy_string(&name,"T1ubq");

    printf("why why");

/*
    error = get_custom_module(&rbasinsmod, rbasinsname);
    if (error != NO_ERROR)
    {
        print_error(error, "init_dsa", 1, " when getting rbasins module");
    }
    free(rbasinsname);


    // Initializing rbasins module.
    execute_init_func_custom_module(rbasinsmod, rbasinscfg);
    free(rbasinscfg);
*/

    // Calculating initial assignement of torsional angles for the whole system.
    update_all_torsional_angles(residues, chains, nchains, coords);

    printf("when when");

    // Initializing global variables..
    T = 350;
    create_bool_array(&flags, nres);


    create_float_array(&curconf, 4 * natoms);
    create_float_array(&optconf, 4 * natoms);


    printf("who who");

    totalnumbins = (360 / par.binsize) * (360 / par.binsize);
//    printf("using totalnumbins: %i\n",totalnumbins);
//
    create_float_array(&rama_distribution, nres * totalnumbins);
    create_index_array(&rama_index, nres * totalnumbins);

    totalpsiphibins = 20 * 20 * 18 * 18;

    append_string(&currentstr,pardir,"/omega_numeric.par");
    create_float_array(&omega_probs, 20 * 20);
    load_omega_probs_from_file(currentstr,omega_probs);

    append_string(&currentstr,pardir,"/psi_psi_numeric.par");
    create_float_array(&psi_psi, totalpsiphibins);
    load_phi_psi_dist_from_file(currentstr,psi_psi);

    append_string(&currentstr,pardir,"/psi_phi_numeric.par");
    create_float_array(&psi_phi, totalpsiphibins);
    load_phi_psi_dist_from_file(currentstr,psi_phi);

    append_string(&currentstr,pardir,"/phi_phi_numeric.par");
    create_float_array(&phi_phi, totalpsiphibins);
    load_phi_psi_dist_from_file(currentstr,phi_phi);

    append_string(&currentstr,pardir,"/phi_psi_long_numeric.par");
    create_float_array(&phi_psi_long, totalpsiphibins);
    load_phi_psi_dist_from_file(currentstr,phi_psi_long);


    append_string(&currentstr,ramadir,"/");
    append_string(&currentstr,currentstr,name);
    append_string(&currentstr,currentstr,".rama_index");

    //load_rama_map_from_file("input/T1ubq.rama_index",rama_index,NULL,TRUE);
    load_rama_map_from_file(currentstr,rama_index,NULL,TRUE);

    append_string(&currentstr,ramadir,"/");
    append_string(&currentstr,currentstr,name);
    append_string(&currentstr,currentstr,".rama_map");
    load_rama_map_from_file(currentstr,NULL,rama_distribution,FALSE);
    //load_rama_map_from_file("input/T1ubq.rama_map",NULL,rama_distribution,FALSE);

    free(currentstr);
    free(pardir);
    free(ramadir);
    free(name);

    //clear_float_array(rama_distribution, nres * totalnumbins);
    //clear_index_array(rama_index, nres * totalnumbins);
    //clear_float_array(omega_probs, 20 * 20);

    //clear_float_array(psi_phi, totalpsiphibins);
    //clear_index_array(psi_psi, totalpsiphibins);
    //clear_index_array(phi_psi_long, totalpsiphibins);

    if(loadseed){
        IntValue tmpseed=0;
        FILE* SeedFile=fopen(seedfile,"r");
        if(SeedFile==NULL) loadseed=FALSE;
        else{    
        if(strcmp(seedfile,"/dev/urandom") == 0  || strcmp(seedfile,"/dev/random") == 0) {
            fread(&tmpseed,sizeof(IntValue),1,SeedFile);
        }
        else fscanf(SeedFile,"%i",&tmpseed);
        printf("        Using seed %i from file %s.\n",tmpseed,seedfile); 
        fclose(SeedFile);
        free(seedfile);
        set_seed(tmpseed);
        }
    }
    if(!loadseed) set_seed(time(0)); 

    if(par.randominit){ 
	if(!restosample){
	generate_rand_conf(); 
   	// printf("\n        Generating inital conformation...");
        // save_opt_conf(); // Current conformation is saved as optimal.
         //printf("done\n");}
	}
	else{
	printf("Generating selected residues");
	generate_rand_selected();
	}
	printf("\n        Generating inital conformation...");
         save_opt_conf(); // Current conformation is saved as optimal.
         printf("done\n");}   
    else {
    	printf("zone1");
    	update_all_coords(coords, distmask, atoms, natoms, residues, nres, chains, nchains, TRUE);
        calculate_all_distances(distances, coords, natoms);
        set_distmask_intervals_product(1,distmask,atoms,natoms,0,natoms-1,0,natoms-1);
        save_opt_conf();
    }

	printf("zone2");

    // Initializing variables.
    curenergy=0;
    for (i=0;i<nenergymods;i++){
        execute_pre_func_energy_module(energymods[i]);
        execute_calc_energy_func_energy_module(energymods[i], &curenergies[i]);
        execute_post_func_energy_module(energymods[i]);
        curenergy+=energycoeffs[i]*curenergies[i];
    }
    printf("            (Initial Energy: %f)\n",curenergy);
    printf("                    ");
    for(i=0;i<nenergymods;i++) printf("  %i) %f*%4.2f  ",i,curenergies[i],energycoeffs[i]);
    printf("\n");


    optenergy = curenergy;
 
    printf("        Total Initialization Time: %i\n",time(NULL)-starttime);
    printf("    Done.\n");

    for (i=0;i<nenergymods;i++) free(energycfgs[i]);
    free(energycfgs);
    energycfgs=NULL;

    

    return error;
}

ErrorCode finish_dsa(void)
{
    IndexValue i;
    delete_float_array(&omega_probs, 20 * 20);
    delete_float_array(&psi_phi, totalpsiphibins);
    delete_float_array(&psi_psi, totalpsiphibins);
    delete_float_array(&phi_phi, totalpsiphibins);
    delete_float_array(&phi_psi_long, totalpsiphibins);

    delete_index_array(&rama_index, nres * totalnumbins);
    delete_float_array(&rama_distribution, nres * totalnumbins);

    delete_float_array(&curconf, 4 * natoms);
    delete_float_array(&optconf, 4 * natoms);

    delete_bool_array(&flags, nres);


//    execute_finish_func_custom_module(rbasinsmod);
    for(i=0;i<nenergymods;i++) execute_finish_func_energy_module(energymods[i]);
    delete_index_array(&energymods,maxnumenergymods);
    delete_float_array(&energycoeffs,maxnumenergymods);
    delete_float_array(&curenergies,maxnumenergymods);
    free(trajectoryname);
    return NO_ERROR;


    free(restomovearray);
}

ErrorCode pre_dsa_step(void)
{
    if (roundnum>0) {  //set local follow up options
        par=followuppar;
        par.countermoves=TRUE;
    }
    par.logupdtempd = log(1 + par.updtempd);
    temp=par.deftemp0;

    //aashish:also set the shared temp
    T=temp;

    energysdevsmooth = par.defenergydev0;
    return NO_ERROR;
}

ErrorCode post_dsa_step(void)
{
    printf("    DSA Run Finished.\n");
    printf("    Total Function Evaluations: %i\n",dsastep );
    printf("    Accepted transitions: %i\n",nincsteps+ndecsteps);
    printf("        Increasing transitions: %i\n",nincsteps);
    printf("        Decreasing transitions: %i\n",ndecsteps);
    printf("    Rejected transitions: %i\n",nrejsteps);
    printf("    Final Energy: %f\n",optenergy);
    printf("    Final Temp: %f\n",temp);
    printf("    Total Running Time: %i\n",time(NULL)-starttime);

    //preparing for a possible second round
    reset_status_parameters();
    roundnum++;
    return NO_ERROR;
}

ErrorCode write_current_conf(FILE* path)
{
    IndexValue a;
    char aa[5];
    char at[5];
    char chain;
    for(a=0;a<natoms;a++)
    {
        chain = 'A' + (char) atoms[a].chain;
        get_aminoacid_name(residues[atoms[a].res].type,aa);
        get_atom_name(atoms[a].type,at);
//         1         2         3         4         5         6         7         8
//	 12345678901234567890123456789012345678901234567890123456789012345678901234567890
//	 MODEL        1
//	 ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N
//	 ATOM      2  CA  ALA A   1      11.639   6.071  -5.147  1.00  0.00           C
        fprintf(path, "ATOM  %5i  %-3s %s %c %3i     %7.3lf %7.3lf %7.3lf  1.00  0.00\n"
	        ,a+1,at,aa,chain,atoms[a].res+1,*atoms[a].x,*atoms[a].y,*atoms[a].z);
    }
}

ErrorCode run_dsa_step(void)
{
    starttime=time(NULL);
    FILE* trjfile;

    if(writetrajectory)
        trjfile = fopen(trajectoryname,"w");
    while (!finished){

    IndexValue res,i,chn;
    BoolValue newconf;
    FloatValue newenergy, denergy;
    FloatValue energymean, energysdev;
  

    res = 0;
    chn = 0;
    newconf = generate_new_conf(&res, &chn);

    if (newconf)
    {
        dsastep++;
        
        newenergy=0;
        for(i=0;i<nenergymods;i++) {
            execute_calc_energy_func_energy_module(energymods[i],&curenergies[i]);
            newenergy+=energycoeffs[i]*curenergies[i];
        }
        denergy = newenergy - curenergy;
        //printf("New Energy %f, CurEnergy %f, OptEnergy %f, Energy change: %f\n", newenergy, curenergy, optenergy, denergy);
        
	if ((denergy <= 0) || (get_frandom() < exp( -denergy / temp)))
        {

          //aashish
	 // printf("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA \n");

            if (denergy < 0) ndecsteps++;
            else nincsteps++;
            

	    //aashish
	    //print_molecule_data(atoms,natoms,residues,nres,chains,nchains);
   
                 if(writetrajectory && saveminimal)
		{
		    fprintf(trjfile, "MODEL %i\n",dsastep);
		    write_current_conf(trjfile);
		    fprintf(trjfile, "ENDMDL\n");
		}
         

	    curenergy = newenergy;
	    
	    //printf("Current Energy is %f \n", curenergy);

	    
	    if (curenergy < optenergy)
            {
                // A new optimal conformation was found.
                optenergy = curenergy;
                
		save_opt_conf();
               if(writetrajectory && saveoptimal)
		{
		    fprintf(trjfile, "MODEL %i\n",dsastep);
		    write_current_conf(trjfile);
		    fprintf(trjfile, "ENDMDL\n");
		}
                printf("        New optimal energy found (T=%f, dsastep=%i, time=%i): %f\n", temp, dsastep,time(NULL)-starttime,optenergy);
                printf("        ");
                for(i=0;i<nenergymods;i++) printf("  %i) %f*%4.2f  ",i,curenergies[i],energycoeffs[i]);
                printf("\n");
            }
        }
        else
        {
	//aashish
	//printf("RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR \n");

	    nrejsteps++;
            restore_prev_conf(res, chn);
        
 	update_all_torsional_angles(residues, chains, nchains, coords);	
	
	}
        
        // Updating energy counters with the energy of the current state.
        nsolutions++;
        energysum += curenergy;
        energysumsq += curenergy * curenergy;
    }

    if (dsastep % par.updtempint == 0)
    {
        // Calculating energy deviation at current temperature.
        energymean = 0.0; 
        energysdev = 0.0;
        calculate_mean_and_sdev(&energymean, &energysdev);

        // Update temperature.
        prevtemp = temp;

        //aashish: also update the shared temp
	T=temp;

        if (0 < energysdevsmooth) temp /= 1.0 + (temp * par.logupdtempd) / (3 * energysdevsmooth);

        // Calculating smoothed deviation using exponential formula (see "The Annealing Algorithm", RHJM Otten, LPPP van Ginneken,
        // page 174, implementation of "smooth" function).
        energysdevsmooth = (1 - par.energyrmsdsmooth) * energysdev + par.energyrmsdsmooth * energysdevsmooth * (temp / prevtemp);

        if (energysdevsmooth <= 0) convergence = TRUE; // This situation should be reported as an error, becase the standard deviation should never get negative.

        if (par.checktempconv)
        {
            // Checking convergence on temperature. The annealing is outside the linear area.
            if (0 < temp) convergence = energysdevsmooth * energysdevsmooth <= par.convepsilon * temp;
            else convergence = TRUE; // When we reach 0 temperature (if this ever happens), we should be at the global minimum (yeah, right).
        }

        // Clean counters.
        nsolutions = 0;
        energysum = 0;
        energysumsq = 0;
    }


    //printf("X: %f Y:%f Z%f for atom 50\n",coords[xcoord_index(50)],coords[ycoord_index(50)],coords[zcoord_index(50)]);
    finished=(convergence || dsastep == par.maxnumsteps);
    if (finished) restore_opt_conf();
    }

/*
    FloatValue energy;
    execute_calc_energy_func_energy_module(energymod, &energy);
    printf("    DOPE energy: %f\n", energy);
    T -= 10.0;

    printf("    Flags array: ");
    IndexValue i;
    for (i = 0; i < min(62, nres); i++)
        printf("%d", flags[i]);
    printf("\n");

    printf("    Just a random number from the Mersenne RNG: %f\n", get_frandom());

    execute_run_funcf_custom_module(rbasinsmod, NULL, 0);

*/


    if(writetrajectory)
        fclose(trjfile);

//    printf(" \n Trajectory saved in  %s \n", trajectoryname);



    return NO_ERROR;
}
