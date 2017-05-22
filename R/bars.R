# Copyright 2014 Patrick O. Perry
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


singlebars <- function(term)
{
    if (is.name(term) || !is.language(term))
        return(term)
    if (length(term) == 2) {
        term[[2]] <- singlebars(term[[2]])
        return(term)
    }
    stopifnot(length(term) >= 3)
    if (is.call(term) && term[[1]] == as.name("||"))
        term[[1]] <- as.name("|")
    for (j in 2:length(term)) term[[j]] <- singlebars(term[[j]])
    term
}


order_bars_hierarchy <- function(bars)
{
    bars_char <- Map(function(bar) as.character(bar)[3], bars)
    # remove parenthesis
    bars_char <- Map(gsub, bars_char, pattern = c("[\\(, \\)]"),
                     replacement = "")
    bars_groups <- Map(function(x) unlist(strsplit(x, split = ":")), bars_char)

    sorted_group_counts <- sort(table(unlist(bars_groups)), decreasing = TRUE)
    if (!all(sorted_group_counts == rev(seq_along(sorted_group_counts))))
        stop("random effects formula misspecified ")

    bars[order(sapply(bars_groups, length))]
}
