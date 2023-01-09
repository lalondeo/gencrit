#include "gtools.h"
#include <stdlib.h>
#include <math.h>
#include "assert.h"

#define MAX_NUMBER_OF_CLIQUES 50
#define MAX_COLORINGS 200
#define MAX_COLORING_ANTICLIQUES 200
#define MAX_WEIGHT_COLORING_ANTICLIQUES 2
#define NUMBER_OF_CONTRACTIONS 4
#define MAX_CLIQUE_SIZE 10
#define MAX_CHI TARGET_CHI+2


void compute_degrees(const graph * g, int n, int * degs)
{
	for(int i = 0; i < n; i++)
		degs[i] = __builtin_popcount(g[i]);
}	

// Turns the binary number a_1 a_2 ... a_{i-1} a_i a_{i+1} ... into a_1 a_2 ... a_{i-1} a_{i+1} ...
setword pop_bit(setword x, int i)
{
	int y = (x << i) >> i;
	return ((y & ~bit[i]) << 1) | (x ^ y);
}
setword expand_bit(setword x, int i)
{
	int y = (x << i) >> i;
	return (y >> 1) | (x ^ y);
}



#define E(g,i,j) g[i] & bit[j]

int find_largest_clique(const graph * g, int n, setword * _clique, int lower_bound)
{
	int vertices[MAXN];
	int i_vertices = 0;
	int current_vertex = 0;
	setword clique = 0;
	int clique_size = 0;
	int omega = lower_bound - 1;

	int degs[MAXN];
	compute_degrees(g, n, degs);
	
	while(1)
	{

		if(current_vertex >= n)
		{
			if(i_vertices == 0) break;
			current_vertex = vertices[--i_vertices];
			
	
			clique_size--;
			clique ^= bit[current_vertex];
			
		}
		
		else if(degs[current_vertex] >= omega-1 && !(clique & ~(g[current_vertex])))
		{
			clique_size++;
			vertices[i_vertices++] = current_vertex;
			clique ^= bit[current_vertex];	
			if(clique_size > omega)
			{
				if(_clique != NULL)
					*_clique = clique;
					
				omega = clique_size;
			}
	
		}

		current_vertex++; 
		
	}
	
	#ifdef DEBUG
	clique = *_clique;
	assert(__builtin_popcount(clique) == omega);
	for(int i = 0; i < n; i++)
		if(clique & bit[i])
			for(int j = i+1; j < n; j++)
				if(clique & bit[j])
					assert(E(g,i,j));
					
	#endif
	
	return omega;
}


boolean test_clique(graph * g, int n, setword clique)
{
	
	for(int i = 0; i < n; i++)
		if((clique & bit[i]) && ((clique ^ bit[i]) & ~(g[i])))
			return FALSE;
			
	return TRUE;
	
}

// Assumes that clique is not the whole graph
// Takes as arguments a graph g, the size n, the number of colors in the coloring k, a clique and a callback function for when a coloring is found
// If said callback function returns FALSE, the function aborts the search
// Returns whether the search was completed (TRUE) or interrupted by the callback function (FALSE)
boolean list_colorings(const graph* g, int n, int k, setword clique, boolean (*log_coloring)(const setword*))
{
	#ifdef DEBUG
	assert(test_clique(g,n,clique));
	#endif
	if(__builtin_popcount(clique) > k) return TRUE;
	

	int best_score, score, i;
	int current_vertex, vertex, current_color;
	int current_bit;
	int colors[MAXN];
	int vertices[MAXN];
	int i_vertices = 0;
	int target = n - __builtin_popcount(clique) - 1; // Corresponds to the value of i_vertices for which we can conclude that every vertex was assigned a color
	
	/*** Setup of the initial state of coloring ***/
	// coloring[i] corresponds to the vertices that were assigned color i so far (as a bit set)
	// This code assigns color i to the i'th vertex in the clique and sets every other vertex as uncolored
	setword coloring[MAX_CHI];
	boolean is_free[MAXN];

	int c = 0;
	for(int i = 0; i < n; i++)
	{
		current_bit = bit[i];
		if(clique & current_bit)
		{
			coloring[c++] = current_bit;
			is_free[i] = FALSE;
		}
		
		else
			is_free[i] = TRUE;
	}
	for(; c <= k; c++)
		coloring[c] = 0;
	
	if(__builtin_popcount(clique) == n)
	{
		(*log_coloring)(coloring);
		return TRUE;
	}


	/*** Branching phase ***/
	choose_new_vertex: 
		// For when we need to pick a value for current_vertex
		// We look for the free vertex with the most colors in its neighborhood, with ties broken according to precedence
		// In a real DSATUR branch-and-bound algorithm, ties would be broken according to the largest number of free vertices in a neighborhood, but this is too expensive for the graph sizes we're interested in
		if(i_vertices <= 3)
		{
			best_score = k+1; // k - best DSATUR found so far
			for(vertex = 0; vertex < n; vertex++)
				if(is_free[vertex])
				{
					score = k;
					for(i = 0; i < k; i++)
						score -= ((g[vertex] & coloring[i]) != 0); // 1 iff there is some vertex in the neighborhood of the vertex with color i
					
					if(score < best_score)
					{
						if(score == 0) // A vertex has all k colors in its neighborhood, and therefore can't be colored
							goto backtrack;
							
						current_vertex = vertex;
						best_score = score;
					}
				}
		}
		else
		{
			for(vertex = 0; vertex < n; vertex++)
				if(is_free[vertex])
				{
					current_vertex = vertex;
					break;
				}
		}
		
		is_free[current_vertex] = FALSE;
		current_color = 0;
		goto try_color;
		
	backtrack:
		if(i_vertices == 0) // No lower level, meaning that the search space has been exhausted
			return TRUE;
	
		i_vertices--;
		current_color = colors[i_vertices];
		current_vertex = vertices[i_vertices];
		coloring[current_color++] ^= bit[current_vertex]; // Sets the current vertex as not colored with colors[i_vertices] anymore
	
	try_color:
		
		if(current_color >= k)
		{
			is_free[current_vertex] = TRUE;
			goto backtrack;
		}
		else if((coloring[current_color] & g[current_vertex]) == 0) // Checks if colors[i_vertices] is a legal color to assign to the current vertex; if so, go to the next level
		{	
			if(i_vertices != target) // Not everything has been colored
			{				
				coloring[current_color] ^= bit[current_vertex]; // Set current vertex as having color colors[i_vertices]
				colors[i_vertices] = current_color;
				vertices[i_vertices++] = current_vertex;
				goto choose_new_vertex;
			}
			coloring[current_color] ^= bit[current_vertex];
			if(!(*log_coloring)(coloring)) 
				return FALSE;
				
			coloring[current_color] ^= bit[current_vertex];	
		}
		current_color++;
		goto try_color;
		
}

boolean callback_has_coloring(const setword* coloring)
{
	return FALSE;
}

boolean has_coloring(const graph* g, int n, char k, setword clique)
{
	int size_clique = __builtin_popcount(clique);
	if(size_clique > k) return FALSE;
	else if(size_clique == n) return TRUE;
	
	return !list_colorings(g, n, k, clique, &callback_has_coloring);
}


void test_coloring(const graph * g, int n, const setword * coloring, int k)
{
	int coloring2[MAXN];
	for(int i = 0; i < n; i++)
	{
		int c = -1;
		for(int j = 0; j < k; j++)
			if((coloring[j] & bit[i]) != 0)
			{
				if(c != -1) {assert(0); }
				c = j;
			}
		if(c == -1) {assert(0);}
		coloring2[i] = c;
	}
	for(int i = 0; i < n; i++)
		for(int j = i+1; j < n; j++)
			if(E(g,i,j) && (coloring2[i] == coloring2[j]))
				assert(0);
				
}



/******** Contracted colorings ********/

setword contracted_colorings[NUMBER_OF_CONTRACTIONS][MAX_COLORINGS][MAX_CHI];
int i_contractions;
int n_contracted_colorings[NUMBER_OF_CONTRACTIONS];
int contraction_order[NUMBER_OF_CONTRACTIONS];


// Assumes that there's an edge between e1 < e2
// Removes vertex e2 and puts edges between e1 and vertices that shared an edge with e2
void contract_edge(const graph * gi, graph * gf, int n, int e1, int e2)
{
	int i;
	for(i = 0; i < e2; i++)
		gf[i] = pop_bit(gi[i], e2);
	
	for(; i < n-1; i++)
		gf[i] = pop_bit(gi[i+1], e2);
		
	gf[e1] |= pop_bit(gi[e2] ^ bit[e1], e2);
	
	for(i = 0; i < n-1; i++)
		if(gf[e1] & bit[i])
			gf[i] |= bit[e1];
}

graph contracted_graph[MAXN];

boolean contracted_colorings_has_all_colorings;
boolean contracted_colorings_callback(const setword * coloring)
{
	if(n_contracted_colorings[i_contractions] == MAX_COLORINGS)
	{
		contracted_colorings_has_all_colorings = FALSE;
		return FALSE;
	}

	for(int k = 0; k < TARGET_CHI-1; k++)
		contracted_colorings[i_contractions][n_contracted_colorings[i_contractions]][k] = coloring[k];
	
	n_contracted_colorings[i_contractions]++;
	return TRUE;
}

boolean generate_contracted_colorings(const graph * g, int n, setword _clique)
{
	setword clique;
	boolean found_clique;
	i_contractions = 0;

	
	#ifdef DEBUG
	int colors[MAXN];
	#endif
	
	for(int i = 1; i < n; i++)
	{
		clique = _clique;
		if(clique & bit[i])
			clique ^= bit[i];
			
		
		for(int j = 0; j < i; j++)
			if(E(g,i,j))
			{
				n_contracted_colorings[i_contractions] = 0;
				contract_edge(g, contracted_graph, n, j, i);
				contracted_colorings_has_all_colorings = TRUE;
				
				list_colorings(contracted_graph, n-1, TARGET_CHI-1, pop_bit(clique, i), &contracted_colorings_callback);
				
				if(n_contracted_colorings[i_contractions] == 0) // Graph is not TARGET_CHI-1-colorable, which is bad
					return FALSE;
					
				else if(contracted_colorings_has_all_colorings) // Bring i back in the story
				{
					for(int l = 0; l < n_contracted_colorings[i_contractions]; l++)
					{
						#ifdef DEBUG
						for(int s = 0; s < n; s++) colors[s] = -1;
						#endif
						
						for(int k = 0; k < TARGET_CHI-1; k++)
						{
							contracted_colorings[i_contractions][l][k] = expand_bit(contracted_colorings[i_contractions][l][k], i);
							if(contracted_colorings[i_contractions][l][k] & bit[j])
								contracted_colorings[i_contractions][l][k] ^= bit[i];
							
							#ifdef DEBUG
							for(int s = 0; s < n; s++)
								if(bit[s] & contracted_colorings[i_contractions][l][k])
								{
									assert(colors[s] == -1);
									colors[s] = k;
								}
							#endif
						}
						
						#ifdef DEBUG
						assert(colors[i] == colors[j]);

						for(int i_ = 1; i_ < n; i_++)
							for(int j_ = 0; j_ < i_; j_++)
								if(E(g,i_,j_))
									assert((colors[i_] == colors[j_]) == ((i == i_) && (j == j_)));
								
						#endif
					}
	
				
					if(++i_contractions == NUMBER_OF_CONTRACTIONS)
						return TRUE;
				}
			}
	
	}
	return TRUE;
}

int cmp_inv_popcount (const void * _a, const void * _b) {
   setword a = (*(setword*)_a);
   setword b = (*(setword*)_b);
   return __builtin_popcount(b) - __builtin_popcount(a);
}

// Computed for the graph corresponding to the first nmax-1 vertices (which must have chromatic number TARGET_CHI-1)
// Corresponds to a list of small independent sets such that it is possible to find a coloring of the whole graph such that, for some given color, the vertices colored with that color are exactly the vertices in the anticlique
// Makes it easy to test for whether the full graph has chromatic number TARGET_CHI-1, by checking whether one of the 
setword _coloring_anticliques[MAX_WEIGHT_COLORING_ANTICLIQUES][MAX_COLORING_ANTICLIQUES];
int _n_coloring_anticliques[MAX_WEIGHT_COLORING_ANTICLIQUES];

setword coloring_anticliques[MAX_COLORING_ANTICLIQUES * MAX_WEIGHT_COLORING_ANTICLIQUES];
int n_coloring_anticliques;


// Called from list_colorings from make_coloring_words when maxn = n + 1
boolean log_colorings_chrom(const setword * coloring)
{
	for(int i = 0; i < TARGET_CHI-1; i++)
	{
		setword anticlique = coloring[i];
		int weight = __builtin_popcount(anticlique);
		if(weight <= MAX_WEIGHT_COLORING_ANTICLIQUES)
		{
			weight--;
			if(_n_coloring_anticliques[weight] == MAX_COLORING_ANTICLIQUES) { printf("RIP %d\n", weight); continue; }
			
			// Makes sure that the anticlique is not already known
			for(int j = 0; j < _n_coloring_anticliques[weight]; j++)
			{
				if(anticlique == _coloring_anticliques[weight][j])
				
					goto continue_loop;
			}	
			 _coloring_anticliques[weight][_n_coloring_anticliques[weight]++] = anticlique;
		}
		
		continue_loop: ;
	}
	return TRUE;
		
}

void make_coloring_anticliques(const graph * g, int n, setword clique)
{
	n_coloring_anticliques = 0;
	for(int i = 0; i < MAX_WEIGHT_COLORING_ANTICLIQUES; i++)
		_n_coloring_anticliques[i] = 0;
	
	list_colorings(g, n, TARGET_CHI-1, clique, &(log_colorings_chrom));
	for(int w_1 = 0; w_1 < MAX_WEIGHT_COLORING_ANTICLIQUES; w_1++)
	{
		for(int i_1 = 0; i_1 < _n_coloring_anticliques[w_1]; i_1++)
		{
			
			// Looks to see if, for some already known anticlique a_, it holds that !(a_ & ~a), where a = _coloring_anticliques[w_1][i_1] is the current clique, and if so, a is rejected
			// This is because this means that for any word b, !(b & a_) -> !(b & a), so that a is useless in proving that a graph has chromatic number TARGET_CHI-1

			for(int w_2 = 0; w_2 < w_1; w_2++)
				for(int i_2 = 0; i_2 < _n_coloring_anticliques[w_2]; i_2++)
					if(!(_coloring_anticliques[w_2][i_2] & ~_coloring_anticliques[w_1][i_1])) 
						goto next;
						
			coloring_anticliques[n_coloring_anticliques++] = _coloring_anticliques[w_1][i_1];
			
			next: ;
		}
	}
			
}

#define max(x,y) (x<y ? y : x)

boolean preprune_crit(graph * g, int n, int maxn)
{
	if(n == maxn)
	{
		//return 1;
		setword lastrow = g[n-1];
		for(int i = 0; i < n_coloring_anticliques; i++)
			if(!(lastrow & coloring_anticliques[i]))
				return TRUE;	
			
	
		for(int i = 0; i < i_contractions; i++)
		{
			for(int j = 0; j < n_contracted_colorings[i]; j++)
				for(int color = 0; color < TARGET_CHI - 1; color++)
					if(!(lastrow & contracted_colorings[i][j][color]))
						goto continue_loop;
						
			return TRUE;	
			continue_loop: ;
		}
	}
	return FALSE;
}

setword base_clique;
int base_chromatic_number;
int base_omega;

boolean prune_crit(graph *g, int n, int maxn)
{

	if(n == maxn)
	{
		
		setword clique;
		setword _clique;
		int omega = find_largest_clique(g, n, &clique, 1);
		int _omega;
		if(omega >= TARGET_CHI && n > omega)
			return TRUE;
			
		if(has_coloring(g, n, TARGET_CHI - 1, clique))
			return TRUE;
		
		graph g2[MAXN];
		for(int i = 0; i < n; i++)
			g2[i] = g[i];
			
		for(int i = n-1; i >= 0; i--)
			for(int j=i-1; j >= 0; j--)
				if(E(g,i,j))
				{
					g2[i] ^= bit[j];
					g2[j] ^= bit[i];
					// Cook up a clique for g2
					_omega = find_largest_clique(g2, n, &_clique, omega-2);
					if(!has_coloring(g2, n, TARGET_CHI-1, _clique))
						return TRUE;
					
					g2[i] ^= bit[j];
					g2[j] ^= bit[i];
				}
	
		return FALSE;
			
	}
	
	else if(n+1 == maxn)
	{
		
		setword clique;
		int omega;
		if(base_omega == TARGET_CHI-1)
		{
			clique = base_clique;
			omega = base_omega;
		}
		else
		{
			omega = find_largest_clique(g, n, &clique, base_omega);
			if(omega >= TARGET_CHI)
				return TRUE;
		}
	
		
		int chi = max(omega, base_chromatic_number);
		if(has_coloring(g, n, chi, clique))
		{
			if(chi + (maxn - n) < TARGET_CHI)  return TRUE;
		}
		else if(++chi == TARGET_CHI)
			return TRUE; 
		


		if(!generate_contracted_colorings(g, n, clique)) return TRUE;
		make_coloring_anticliques(g, n, clique);	
		return FALSE;
		
	}
	else if(n+2 == maxn)
	{
		base_omega = find_largest_clique(g, n, &base_clique, 1);
		if(base_omega >= TARGET_CHI) return TRUE;
		base_chromatic_number = base_omega;
		while(!has_coloring(g, n, base_chromatic_number, base_clique))
		{
			if(base_chromatic_number == TARGET_CHI-1) return TRUE;
			base_chromatic_number++;
		}
	
		if(base_chromatic_number + (maxn - n) < TARGET_CHI)
			return TRUE;	

		return FALSE;
	
	}
	else
		return FALSE;
	

}

