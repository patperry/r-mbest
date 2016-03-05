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
#include "rmbest.h"


SEXP rdglm_index_loop(SEXP sgroup, SEXP subset)
{
	SEXP subset_g;
	int j, k, n = LENGTH(sgroup);
	int g, ngroups = LENGTH(subset);
	const int *group = INTEGER(sgroup);
	int *group_size = (void *)R_alloc(ngroups, sizeof(*group_size));

	for (g = 0; g < ngroups; g++) {
		group_size[g] = 0;
	}

	for(k = 0; k < n; k++) {
		g = group[k] - 1; //convert from R to C indexing
		j = group_size[g];
		subset_g = VECTOR_ELT(subset, g);
		INTEGER(subset_g)[j] = k + 1; // output R indexing
		group_size[g] = j + 1;
	}

	return R_NilValue;
}
