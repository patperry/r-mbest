/*
 * Copyright 2016 Patrick O. Perry.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <Rdefines.h>
#include <string.h>
#include "rmbest.h"


SEXP group_subsets(SEXP sgroup, SEXP sngroups)
{
	SEXP subset, subset_g;
	const int *group;
	R_len_t *group_size;
	R_len_t i, j, nobs;
	int g, sg, ngroups;

	/* validate input */
	if (TYPEOF(sgroup) != INTSXP) {
		error("invalid 'group' argument");
	}
	group = INTEGER(sgroup);
	nobs = XLENGTH(sgroup);

	ngroups = asInteger(sngroups);
	if (ngroups == NA_INTEGER || ngroups < 0) {
		error("invalid 'ngroups' argument");
	}

	/* determine group sizes */
	group_size = (void *)R_alloc(ngroups, sizeof(*group_size));
	memset(group_size, 0, ngroups * sizeof(*group_size));

	for (i = 0; i < nobs; i++) {
		sg = group[i];
		if (sg == NA_INTEGER || sg <= 0) {
			/* ignore observations with invalid group indices */
			continue;
		}
		g = sg - 1;
		group_size[g] += 1;
	}

	/* allocate group index vectors, reset group sizes to 0 */
	PROTECT(subset = allocVector(VECSXP, ngroups));
	for (g = 0; g < ngroups; g++) {
		/* use REALSXP instead of INTSXP to allow indices > 2^31 */
		SET_VECTOR_ELT(subset, g, allocVector(REALSXP, group_size[g]));
		group_size[g] = 0;
	}

	/* copy group index vectors */
	for (i = 0; i < nobs; i++) {
		sg = group[i];
		if (sg == NA_INTEGER || sg <= 0) {
			continue;
		}
		g = sg - 1;
		subset_g = VECTOR_ELT(subset, g);

		j = group_size[g];
		REAL(subset_g)[j] = (double)(i + 1);
		group_size[g] = j + 1;
	}

	UNPROTECT(1);
	return subset;
}
