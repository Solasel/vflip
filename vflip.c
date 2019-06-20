/* TODO */
/* Incorporate the value of the next game. This will
   require some info about expected winrates, etc. */
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>

#define LOG 0

#define ALL_VALS (15U)

/* Returns a bitmap identical to d but with the xth bit set. */
#define SET_BIT(d, x) ((d) | (1U << (x)))

/* Same as above, but with the xth bit reset. */
#define RESET_BIT(d, x) ((d) & ~(1U << (x)))

/* Returns nonzero if the xth bit is set, and zero if reset. */
#define GET_BIT(d, x) ((d) & (1U << (x)))

/* Returns a true value if a domain is dangerous to click on.
   IE: can this tile be a voltorb? */
#define UNSAFE(x) ((x) & 1U)

/* Returns the row/column corresponding
   to the given variable x. */
#define ROW(x) ((int) (x) / 5)
#define COL(x) ((int) (x) % 5)

/* Returns the variable corresponding
   to row r and column c. */
#define VAR(r, c) (5 * (r) + (c))

/* Gets the constraint corresponding to n variables which
   sum to s. See below for more info. */
#define GET_PART_CONS(n, s) (part_cons[6 * (s) + (n)])

/* Returns the ith solution's vth variable as stored in s. */
#define SOLN_ASGMT(s, i, v) (s[25 * (i) + (v)])

/* Returns the probability of variable v having value x as stored in p. */
#define PROB(p, v, x) (p[4 * (v) + (x)])

/* Maximum depth of tree search. */
#define MAX_DEPTH 4

/* Macro for the evaluation function so we can easily change it. */
#define EVAL(b, m) score(b, m)

/* A map to the number of possible values a domain can take on. */
const int cardinality[16] = {
	0, 1, 1, 2,
	1, 2, 2, 3,
	1, 2, 2, 3,
	2, 3, 3, 4
};

/* Bitmaps corresponding to valid assignments of
   vars given sum and num of vars. The map
   corresponding to n vars summing to s is
   at part_cons[6 * s + n], given by the macro above. */
const uint32_t part_cons[96] = {
	0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
	0x0, 0x2, 0x0, 0x0, 0x0, 0x0,
	0x0, 0x4, 0x2, 0x0, 0x0, 0x0,
	0x0, 0x8, 0x6, 0x2, 0x0, 0x0,
	0x0, 0x0, 0xE, 0x6, 0x2, 0x0,
	0x0, 0x0, 0xC, 0xE, 0x6, 0x2,
	0x0, 0x0, 0x8, 0xE, 0xE, 0x6,
	0x0, 0x0, 0x0, 0xE, 0xE, 0xE,
	0x0, 0x0, 0x0, 0xC, 0xE, 0xE,
	0x0, 0x0, 0x0, 0x8, 0xE, 0xE,
	0x0, 0x0, 0x0, 0x0, 0xE, 0xE,
	0x0, 0x0, 0x0, 0x0, 0xC, 0xE,
	0x0, 0x0, 0x0, 0x0, 0x8, 0xE,
	0x0, 0x0, 0x0, 0x0, 0x0, 0xE,
	0x0, 0x0, 0x0, 0x0, 0x0, 0xC,
	0x0, 0x0, 0x0, 0x0, 0x0, 0x8
};

/* Linked list structure holding a solution to a given board. */
struct soln_elem {

	struct soln_elem *next;
	int assigns[25];
};

/* Models an instance of a game board. Constraints are stored in
   the first four arrays, but our knowledge about the assignments of
   variables is different. Early on in the exploration process, our
   knowledge is encoded in domains. After enumeration, this array
   becomes deprecated, and instead our information is stored in
   solns, which is a flattened array that contains all possible
   solutions to the constraints given what we know. As we gain more
   info, rather than re-alloc the array, we use mask to indicate
   which solutions are no longer valid. */
struct board {

	/* Row/column constraints given by the board. */
	int row_sums[5];
	int row_volts[5];
	int col_sums[5];
	int col_volts[5];

	/* Represents the assignment state of a variable.
	   -1 for not assigned or the value of the variable.
	   Early on, this array is used for arc consistency/domain
	   checking. Later, it is only maintained for printing. */
	int assigns[25];

	/* A bitmap corresponding to the
	   remaining domains of each variable.
	   Specifically, the nth bit is set if n
	   is a valid value for the variable.
	   No longer valid after enumerating boards. */
	uint32_t domains[25];

	/* How many variables in the row/column
	   have been assigned already? */
	int row_asgmt[5];
	int col_asgmt[5];

	/* Has the partition constraint on a
	   row/column been enforced or not? */
	bool row_cons[5];
	bool col_cons[5];

	/* Has the row/column arc of this variable
	   been enforced? */
	bool row_arc[25];
	bool col_arc[25];

	/* An array to hold all the possible solutions
	   to the board. The array gets populated once,
	   and then instead of removing elements from it,
	   we just set mask[i] to true to invalidate it. */
	int *solns, solns_len;
	bool *mask;

	/* An array holding the probabilities of each
	   value occuring in a certain variable. */
	double probs[100];
};

/* Algorithmic functions. */
void init_board(struct board *b, int *rs, int *rv, int *cs, int *cv);
void solve_board(struct board *b);
bool enf_part(struct board *b);
bool enf_arcs(struct board *b);
bool pfree_doms(struct board *b);
bool enum_boards(struct board *b);
void populate_probs(struct board *b, bool *mask, double *probs);
bool pfree_mask(struct board *b);
bool p_best_act(struct board *b);
uint32_t row_part_cons(struct board *b, int x);
uint32_t col_part_cons(struct board *b, int x);
bool change_var_dom(struct board *b, int var, uint32_t dom);
bool check_row_arc(struct board *b, int var, int val);
bool check_col_arc(struct board *b, int var, int val);
bool r_enum_boards(struct board *b, struct soln_elem **last);
void assign_var_mask(struct board *b, int var, int val, bool *mask);
bool is_solved(struct board *b, bool *mask);
double tree_search(struct board *b, bool *m_in, int var, int val, int depth);
void assign_var_doms(struct board *b, int var, int val);
void incorp_dom(struct board *b, int var, uint32_t *c);
bool ll_append(int *s, struct soln_elem **last);
int score(struct board *b, bool *mask);

/* Testing functions. */
void print_cons(struct board *b);
void print_assigns(struct board *b);
void print_doms(struct board *b);
void print_solns(struct board *b);
void print_probs(double *probs);

/* Main. */
int main()
{
	struct board _game;
	struct board *game = &_game;
	
	/* TODO */
	/* Input. */

#if 0
	/* 3-rooks. */
	int rs[] = {2, 2, 2, 0, 0};
	int rv[] = {4, 4, 4, 5, 5};
	int cs[] = {2, 2, 2, 0, 0};
	int cv[] = {4, 4, 4, 5, 5};
#endif
	int rs[] = {8, 5, 2, 5, 5};
	int rv[] = {1, 1, 3, 2, 3};
	int cs[] = {4, 9, 5, 2, 5};
	int cv[] = {2, 1, 2, 3, 2};
	init_board(game, rs, rv, cs, cv);
	solve_board(game);

	return 0;
}

/* Initializes a board with the given row and column constraints. */
void init_board(struct board *b, int *rs, int *rv, int *cs, int *cv)
{
	int i;
	uint32_t *dom = b->domains;
	int *asgmt = b->assigns;

	/* Copy over all the constraints. */
	memcpy(b->row_sums, rs, 5 * sizeof(int));
	memcpy(b->row_volts, rv, 5 * sizeof(int));
	memcpy(b->col_sums, cs, 5 * sizeof(int));
	memcpy(b->col_volts, cv, 5 * sizeof(int));

	/* Domains of all variables
	   unconstrained and unassigned. */
	for (i = 0; i < 25; i++) {
		dom[i] = ALL_VALS;
		asgmt[i] = -1;
	}

	/* No variables in any row/col have been assigned... */
	memset(b->row_asgmt, 0, 5 * sizeof(int));
	memset(b->col_asgmt, 0, 5 * sizeof(int));

	/* No partition constraints have been enforced... */
	memset(b->row_cons, 0, 5 * sizeof(bool));
	memset(b->col_cons, 0, 5 * sizeof(bool));

	/* And no vars have been arc enforced. */
	memset(b->row_arc, 0, 25 * sizeof(bool));
	memset(b->col_arc, 0, 25 * sizeof(bool));

	/* No solutions computed yet. */
	b->solns_len = 0;
	b->solns = NULL;
	b->mask = NULL;
}

/* Solves a board that has been initialized.

   1. Enumerate all possible boards that satisfy the information
	we have so far, using backtracking search with arc consistency.

   - Now, instead of dealing with constraints, we will deal with possible boards.

   2. Prompt the user to reveal all tiles that are now 100% safe, and remove all
	boards that contradict the new information. Repeat until all remaining
	variables are dangerous.

   3. Now, use MAX_DEPTH game tree analysis to find the best option to continue.
	Specifically, voltorb flip can be modeled as an expectimax game tree with
	two agents - the player, and the randomizer.

	The player makes one of 26 actions on their turn: revealing one of the tiles
	or quitting. If the player chooses to reveal a tile, the randomizer then
	randomly decides the value of that tile from a distribution that can easily
	be recovered from the list of boards.

	Simply do recursive tree search to find the expected value of each action the
	player can take, and choose the best. If the game is over (through winning,
	losing, or quitting), we're done. Otherwise, go back to step 2. */
void solve_board(struct board *b)
{
	/* Reduce using arc consistency, since enum_boards
	   expects an already arc-reduced board. */
	while (enf_part(b) || enf_arcs(b) || pfree_doms(b));

	/* Create a list of all possible boards given the constraints. */
	if (enum_boards(b)) {
		printf("ERROR: Ran out of memory mallocing solution LL elements - %d\n", __LINE__);
		return;
	}


	/* Populate the probabilities so that we can pass
	   this board into pfree_mask(). */
	populate_probs(b, b->mask, b->probs);

	/* Prompt the user to reveal all tiles that are safe,
	   until there are no more safe tiles. Then, do the
	   best action as given by our game theory engine.
	   Repeat until the game is over. */
	while (pfree_mask(b) || p_best_act(b));

	free(b->solns);
	free(b->mask);
}

/* Prompts the user for the true value of all safe tiles on
   the board. Returns true if we make any assignments. This version
   operates on the domains framework, to reduce the amount of
   memory allocated for the solution list. */
bool pfree_doms(struct board *b)
{
	int free_var, input, *asgmts = b->assigns;
	uint32_t dom, *doms = b->domains;

	/* Find the first instance of a safe variable. It must be
	   	1. Unassigned
		2. Safe => no chance of being a 0. */
	for (free_var = 0; free_var < 25; free_var++) {
		if (asgmts[free_var] == -1 && !UNSAFE(doms[free_var]))
			break;
	}

	/* If we couldn't find any variables,
	   we can't make any assignments. */
	if (free_var == 25)
		return false;

	/* Otherwise, prompt the user for the value of the safe variable. */
	print_assigns(b);
	printf("Variable %d is safe. What is its value? It's one of: ", free_var);

	dom = doms[free_var];
	if (GET_BIT(dom, 1))
		printf("1 ");

	if (GET_BIT(dom, 2))
		printf("2 ");

	if (GET_BIT(dom, 3))
		printf("3 ");

	/* Prompts the user for the true value of the tile. */
	printf("\nValue: ");
	while (1) {
		while (scanf("%d", &input) == EOF);

		if (input <= 0 || input > 3 || !GET_BIT(dom, input)) {
			printf("Invalid input, try again.\nValue: ");
			continue;
		}

		break;
	}

	change_var_dom(b, free_var, 1U << input);
	return true;
}

/* Stores all the possible solutions to a board in its solns array. */
bool enum_boards(struct board *b)
{
	int i, num, *solns;
	struct soln_elem ll_sent, *soln_iter, *soln_next, *last;
	bool *mask;

	/* Initialize the linked list sentinel node. */
	ll_sent.next = NULL;
	last = &(ll_sent);

	/* Populates the linked list with all solutions to b. */
	if (r_enum_boards(b, &last))
		return true;

	/* Malloc the (usually large) set of solutions and their mask. */
	num = b->solns_len;
	solns = malloc(25 * num * sizeof(int));
	mask = malloc(num * sizeof(bool));
	if (!solns || !mask) {
		printf("ERROR: Failed to allocate solution/mask buffers: %d.\n", __LINE__);
		return true;
	}

	/* Set the mask to all false. */
	memset(mask, 0, num * sizeof(bool));

	/* Populates the new solns array with the contents of the
	   linked list, and simultaneously frees the linked list. */
	for (i = 0, soln_iter = ll_sent.next; i < num; i++, soln_iter = soln_next) {

		/* Copy the contents of the LL element into the global array... */
		memcpy(solns + 25 * i, soln_iter->assigns, 25 * sizeof(int));

		/* And free the node that we got it from. */
		soln_next = soln_iter->next;
		free(soln_iter);
	}

	if (soln_iter != NULL) {
		printf("ERROR: Iterating through the linked list didn't reach the end: %d.\n", __LINE__);
		exit(0);
	}

	b->mask = mask;
	b->solns = solns;

	return false;
}

/* Populates an array with the probability of each value
   of every variable given a mask for the board b. */
void populate_probs(struct board *b, bool *mask, double *probs)
{
	int i, j, len = b->solns_len, rem = 0, *solns = b->solns;

	/* Clear out the probabilities info... */
	memset(probs, 0, 100 * sizeof(double));

	/* And repopulate them. Since the probability is just
	   num_occurences/total_samples, use the array as a counter
	   for each occurrence, and then divide them all by the
	   remaining number of possible solutions. */
	for (i = 0; i < len; i++) {

		/* If this solution has already been invalidated, ignore it. */
		if (mask[i])
			continue;

		rem++;

		/* For every variable, increment num_occurrences
		   for the value in this solution. */
		for (j = 0; j < 25; j++)
			PROB(probs, j, SOLN_ASGMT(solns, i, j))++;
	}

	/* Finally, divide all the values by the total number of solutions. */
	for (i = 0; i < 100; i++)
		probs[i] /= rem;

}

/* Updates b->assigns to reflect variables that have been
   recently assigned. Assumes b->probs is up to date. */
void update_asgmt(struct board *b)
{
	int i, *asgmt = b->assigns;
	double *probs = b->probs;

	/* For all variables that now have a definite
	   assignment, set b->assigns appropriately. */
	for (i = 0; i < 25; i++) {
		if (PROB(probs, i, 0) == 1.0) {
			asgmt[i] = 0;
			continue;

		} else if (PROB(probs, i, 1) == 1.0) {
			asgmt[i] = 1;
			continue;

		} else if (PROB(probs, i, 2) == 1.0) {
			asgmt[i] = 2;
			continue;

		} else if (PROB(probs, i, 3) == 1.0) {
			asgmt[i] = 3;
			continue;
		}
	}
}

/* Prompts the user to take the best action on the board. Returns true if the game continues. */
bool p_best_act(struct board *b)
{
	int i, j, input, best_a, *asgmt = b->assigns;
	double cand, best = -1.0, *probs = b->probs, q_vals[25];
	uint32_t dom = 0;

	if (is_solved(b, b->mask)) {
		printf("Board is solved!\n");
		print_assigns(b);
		return false;
	}

	/* Find the variable with the largest return for revealing it. */
	for (i = 0; i < 25; i++) {

		q_vals[i] = 0.0;

		/* If the variable is already assigned, we
		   don't gain anything by revealing it. */
		if (asgmt[i] != -1)
			continue;

#if LOG
		printf("Checking variable %d.\n\n", i);
#endif

		/* Otherwise, calculate its expected value. */
		cand  = PROB(probs, i, 1) > 0 ? PROB(probs, i, 1) * tree_search(b, b->mask, i, 1, 1) : 0;
		cand += PROB(probs, i, 2) > 0 ? PROB(probs, i, 2) * tree_search(b, b->mask, i, 2, 1) : 0;
		cand += PROB(probs, i, 3) > 0 ? PROB(probs, i, 3) * tree_search(b, b->mask, i, 3, 1) : 0;

		q_vals[i] = cand;

		/* If it's better, set this as the best action so far. */
		if (cand > best) {
			best = cand;
			best_a = i;
		}
	}

	if (best == -1) {
		printf("Didn't find a next best move at root of search tree.\n");
		exit(0);
	}

#if LOG
	printf("\n\n");
	for (i = 0; i < 5; i++) {
		for (j = 0; j < 5; j++)
			printf("%03.02f ", q_vals[VAR(i, j)]);
		printf("\n");
	}
#endif
	if (score(b, b->mask) > best) {
		print_assigns(b);
		printf("It is not worthwhile to continue. Quitting the game is the best option.\n");
		return false;
	}

	printf("\n\n");
	print_assigns(b);

	printf("The next best action is tile %d. What is its value? It is one of: ", best_a);

	if (PROB(probs, best_a, 0) > 0) {
		dom = SET_BIT(dom, 0);
		printf("0 ");
	}

	if (PROB(probs, best_a, 1) > 0) {
		dom = SET_BIT(dom, 1);
		printf("1 ");
	}

	if (PROB(probs, best_a, 2) > 0) {
		dom = SET_BIT(dom, 2);
		printf("2 ");
	}

	if (PROB(probs, best_a, 3) > 0) {
		dom = SET_BIT(dom, 3);
		printf("3 ");
	}

	/* Prompts the user for the true value of the tile. */
	printf("\nValue: ");
	while (1) {
		while (scanf("%d", &input) == EOF);

		if (input < 0 || input > 3 || !GET_BIT(dom, input)) {
			printf("Invalid input, try again.\nValue: ");
			continue;
		}

		break;
	}

	if (input == 0) {
		printf("zannen.\n");
		return false;
	}

	/* Updates the mask of the board, removing all
	   solns that disagree with input. */
	assign_var_mask(b, best_a, input, b->mask);

	/* Officially assign the value to the variable. */
	asgmt[best_a] = input;

	/* Recalculate the probabilities given the new list of boards. */
	populate_probs(b, b->mask, b->probs);
	update_asgmt(b);

	return true;
}

/* Say we have a board defined by b and m_in. This function either
   returns the expected value of the best action, or, if
   we've hit MAX_DEPTH, the evaluation of the board. */
double tree_search(struct board *b, bool *m_in, int var, int val, int depth)
{
	int _, i, len = b->solns_len;
	double rv = -1.0, cand, probs[100];
	bool *mask = malloc(len * sizeof(bool));
	
	/* Improve OOM handling. */
	if (!mask) {
		printf("Ran out of memory trying to malloc a mask for %d -> %d at depth %d.\n", val, var, depth);
		exit(0);
	}

	/* Update the mask for the result. */
	memcpy(mask, m_in, len);
	assign_var_mask(b, var, val, mask);

#if LOG
	for (_ = 0; _ < depth; _++)
		printf("\t");
	printf("Trying %d -> %d:\n\n", val, var);
#endif

	/* If we're at the maximum depth or have an already
	   solved board, just return the eval function.
	   Otherwise, expand the search to one level deeper. */
	if (depth == MAX_DEPTH || is_solved(b, mask)) {

#if LOG
		for (_ = 0; _ < depth; _++)
			printf("\t");
#endif
		rv = EVAL(b, mask);
#if LOG
		printf("Reached maximum recursion depth! Value of board is %01.02f.\n\n", rv);
#endif


		free(mask);
		return rv;
	}

	/* Calculate the new probability distributions of each variable. */
	populate_probs(b, mask, probs);

	/* Find the variable with the largest return for revealing it. */
	for (i = 0; i < 25; i++) {

		/* If the variable is already assigned, we
		   don't gain anything by revealing it. */
		if (PROB(probs, i, 0) == 1.0 ||
				PROB(probs, i, 1) == 1.0 ||
				PROB(probs, i, 2) == 1.0 ||
				PROB(probs, i, 3) == 1.0)
			continue;

#if LOG
		for (_ = 0; _ < depth; _++)
			printf("\t");
		printf("Checking variable %d.\n\n", i);
#endif

		/* Otherwise, calculate its expected value. */
		cand  = PROB(probs, i, 1) > 0 ? PROB(probs, i, 1) * tree_search(b, mask, i, 1, depth + 1) : 0;
		cand += PROB(probs, i, 2) > 0 ? PROB(probs, i, 2) * tree_search(b, mask, i, 2, depth + 1) : 0;
		cand += PROB(probs, i, 3) > 0 ? PROB(probs, i, 3) * tree_search(b, mask, i, 3, depth + 1) : 0;

#if LOG
		for (_ = 0; _ < depth; _++)
			printf("\t");
		printf("Value of clicking on %d at depth %d is %01.02f.\n\n", i, depth, cand);
#endif

		/* If it's better, set this as the best action so far. */
		if (cand > rv)
			rv = cand;
	}

	free(mask);
	if (rv == -1) {
		printf("Didn't expand any actions for a non-completed board with %d -> %d at depth %d.\n", val, var, depth);
		exit(0);
	}

	return rv;
}

/* Returns the max score obtainable with all known tiles of a
   board b given a mask m. Since this is only called down the
   line in the game theory tree, we can be sure that we know
   at least one tile - our score is nonzero. */
int score(struct board *b, bool *mask)
{
	int i, j, rv = 1, len = b->solns_len, *solns = b->solns;
	uint32_t doms[25];

	memset(doms, 0, 100);

	/* Compute the domains of each variable. */
	for (i = 0; i < len; i++) {
		if (mask[i])
			continue;

		for (j = 0; j < 25; j++)
			doms[j] = SET_BIT(doms[j], SOLN_ASGMT(solns, i, j));
	}

	/* Now iterate over all the domains, and multiply
	   any known vars into the total score. */
	for (i = 0; i < 25; i++) {
		if (doms[i] == 4U)
			rv *= 2;
		else if (doms[i] == 8U)
			rv *= 3;
	}

	return rv;
}

/* Returns whether a given board has been solved - that is,
   whether or not all the hidden 2s and 3s have been found. */
bool is_solved(struct board *b, bool *mask)
{
	int i, j, len = b->solns_len, *solns = b->solns;
	int asgmt[25], soln_asgmt;

	/* Find the first valid board's index. */
	for (i = 0; i < len && mask[i]; i++);

	/* Copy the assignments of the first solution. */
	memcpy(asgmt, solns + 25 * i, 25 * sizeof(int));

	/* Then compare our first possible solution
	   to all the other solutions. */
	for (i++; i < len; i++) {
		if (mask[i])
			continue;

		for (j = 0; j < 25; j++) {
			soln_asgmt = SOLN_ASGMT(solns, i, j);

			/* If the two values differ, the tile is unsure.
			   If then one of the values could be a 2 or a
			   3, we still need to reveal a scoring tile. */
			if (soln_asgmt != asgmt[j] &&
					(soln_asgmt == 2 || soln_asgmt == 3 ||
					 asgmt[j] == 2 || asgmt[j] == 3))
				return false;
		}
	}

	/* Otherwise, we're done. */
	return true;
}

/* Populates the given linked list with all solutions
   that satisfy the board's constraints using backtracking
   search with arc consistency. Returns true if malloc
   fails at some point. Assumes it's passed an already
   arc-reduced board. */
bool r_enum_boards(struct board *b, struct soln_elem **last)
{
	int i, j, asgmts[25], card, min_card = 5, min_asgmt = -1;
	uint32_t doms[25], dom;
	bool valid;

	/* Make copies of the board to backtrack to. */
	memcpy(asgmts, b->assigns, 25 * sizeof(int));
	memcpy(doms, b->domains, 100);

	/* Find the minimum domain unassigned variable to assign. */
	for (i = 0; i < 25; i++) {
		card = cardinality[doms[i]];

		if (asgmts[i] == -1 && card < min_card) {
			min_card = card;
			min_asgmt = i;
		}
	}

	/* If there are no unassigned variables, we have a
	   valid solution to the constraints! Add it to the list. */
	if (min_asgmt == -1) {
		b->solns_len++;
		return ll_append(asgmts, last);
	}

	/* Otherwise, try on each value for the most-constrained
	   variable, and recursively call enum_boards. */
	dom = doms[min_asgmt];
	for (i = 0; i < 4; i++) {

		/* Skip invalid assignments. */
		if (!GET_BIT(dom, i))
			continue;

		/* Assign the variable to i by restricting its domain.
		   This also invalidates arc-consistency. */
		change_var_dom(b, min_asgmt, 1U << i);

		/* Re-enforce arc consistency... */
		while (enf_part(b) || enf_arcs(b));

		/* If any variables have empty domains, this is an invalid assignment. */
		valid = true;
		for (j = 0; j < 25; j++) {
			if (b->domains[j] == 0) {
				printf("Assignment invalid, backtracking.\n");
				valid = false;
			}
		}

		/* Recursively call ourselves, but only if it's a valid assignment. */
		if (valid && r_enum_boards(b, last))
			return true;

		/* And restore the previous state
		   of b for backtracking. */
		memcpy(b->assigns, asgmts, 25 * sizeof(int));
		memcpy(b->domains, doms, 100);
	}

	return false;
}

/* Prompts the user for the true value of all safe tiles on
   the board. Returns true if we make any assignments. Used
   after the solutions to a board are enumerated. */
bool pfree_mask(struct board *b)
{
	int free_var, input, *asgmt = b->assigns;
	double *probs = b->probs;
	uint32_t dom = 0;

	/* Find the first variable to be safe. */
	for (free_var = 0; free_var < 25; free_var++) {

		/* For a variable to be safe, it has to be unassigned
		   and yet have a 0.0% chance to be a 0. */
		if (asgmt[free_var] == -1 && PROB(probs, free_var, 0) == 0.0)
			break;
	}

	/* If there weren't any safe variables, we can't make an assignment. */
	if (free_var == 25)
		return false;

	/* Otherwise, prompt the user for the value of the safe variable. */
	print_assigns(b);
	printf("Variable %d is safe. What is its value? It's one of: ", free_var);

	if (PROB(probs, free_var, 1) > 0) {
		dom = SET_BIT(dom, 1);
		printf("1 ");
	}

	if (PROB(probs, free_var, 2) > 0) {
		dom = SET_BIT(dom, 2);
		printf("2 ");
	}

	if (PROB(probs, free_var, 3) > 0) {
		dom = SET_BIT(dom, 3);
		printf("3 ");
	}

	/* Prompts the user for the true value of the tile. */
	printf("\nValue: ");
	while (1) {
		while (scanf("%d", &input) == EOF);

		if (input <= 0 || input > 3 || !GET_BIT(dom, input)) {
			printf("Invalid input, try again.\nValue: ");
			continue;
		}

		break;
	}

	/* Updates the mask of the board, removing all
	   solns that disagree with input. */
	assign_var_mask(b, free_var, input, b->mask);

	/* Officially assign the value to the variable. */
	asgmt[free_var] = input;

	/* Recalculate the probabilities given the new list of boards. */
	populate_probs(b, b->mask, b->probs);
	update_asgmt(b);

	/* We were able to assign a value, so we win! */
	return true;
}

/* Enforces row and column partition constraints.
   returns true if an assignment was made. */
bool enf_part(struct board *b)
{
	int i, j;
	uint32_t cons;
	bool rv = false;

	/* Iterates through the rows of the board. */
	for (i = 0; i < 5; i++) {

		/* If the row has already been enforced, skip it. */
		if (b->row_cons[i])
			continue;

		b->row_cons[i] = true;

		/* Look up the constraints on the variables. */
		cons = row_part_cons(b, i);

		/* If cons can't change any variables'
		   domains, just go to the next row. */
		if (cons == ALL_VALS)
			continue;

		/* Intersect cons with the domains of all unassigned vars. */
		for (j = 5 * i; j < 5 * i + 5; j++) {

			/* If the variable has been assigned
			   already, ignore it. */
			if (b->assigns[j] >= 0)
				continue;

			/* Otherwise, intersect the new constraints
			   with the variable's current constraints. */
			if (change_var_dom(b, j, b->domains[j] & cons))
				rv = true;
		}
	}

	/* Same as above, but for the columns instead. */
	for (i = 0; i < 5; i++) {

		/* If the col has already been enforced, skip it. */
		if (b->col_cons[i])
			continue;

		b->col_cons[i] = true;

		/* Look up the constraints on the variables. */
		cons = col_part_cons(b, i);

		/* If cons can't change any variables'
		   domains, just go to the next row. */
		if (cons == ALL_VALS)
			continue;

		/* Intersect cons with the domains of all unassigned vars. */
		for (j = i; j < 25; j += 5) {

			/* If the variable has been assigned
			   already, ignore it. */
			if (b->assigns[j] >= 0)
				continue;

			/* Otherwise, intersect the new constraints
			   with the variable's current constraints. */
			if (change_var_dom(b, j, b->domains[j] & cons))
				rv = true;
		}
	}

	/* Return whether or not any variables have been assigned. */
	return rv;
}

/* Enforces arc consistency on all variables in board b.
   Returns as soon as it makes an assignment so enf_part
   can take over and enforce partitioning first. */
bool enf_arcs(struct board *b)
{
	int i, j;
	uint32_t dom;

	/* Iterate over all variables whose arcs haven't been
	   enforced yet and enforce arc consistency. */
	for (i = 0; i < 25; i++) {

		/* If i has already been enforced, skip it. */
		if (b->row_arc[i] && b->col_arc[i])
			continue;

		/* Otherwise, load its domain and
		   try on different valid values. */
		dom = b->domains[i];
		for (j = 0; j < 4; j++) {

			/* If an assignment is invalid, skip it. */
			if (!GET_BIT(dom, j))
				continue;

			/* Otherwise, check that this assignment
			   is arc consistent for rows, cols, or both. */
			if ((!b->row_arc[i] && !check_row_arc(b, i, j)) || (!b->col_arc[i] && !check_col_arc(b, i, j)))
				dom = RESET_BIT(dom, j);
		}

		/* Set the variable's domain to the new, possibly
		   reduced, domain. If we make an assignment because
		   of this restriction, return immediately. */
		if (change_var_dom(b, i, dom))
			return true;
	}

	/* If we get to the end without being able to assign anything, return false. */
	return false;
}

/* Checks the arc from a variable to a row with a certain value.
   Returns true if the value doesn't break arc consistency. */
bool check_row_arc(struct board *b, int var, int val)
{
	int i, row = ROW(var);

	/* GET_BIT(combos[i], j) == 1 iff there is a
	   valid set of assignments to the vars in the
	   row that sums up to i and has j voltorbs.

	   Currently the array is empty, so we
	   construct it in the for loop below. */
	uint32_t combos[16];

	memset(combos + 1, 0, 60);
	combos[0] = 1U;

	if (b->assigns[var] != -1) {
		print_assigns(b);
		printf("Error: checking arc consistency on already-assigned variable %d.\n", var);
		exit(0);
	}

	/* Temporarily assign the variable to the tentative value... */
	b->assigns[var] = val;

	/* Get the possible sums/their voltorbs... */
	for (i = 5 * row; i < 5 * row + 5; i++)
		incorp_dom(b, i, combos);

	/* And reset the variable to not being assigned. */
	b->assigns[var] = -1;

	/* If there's a combination that has the correct row_sum and row_volts, it's valid. */
	return GET_BIT(combos[b->row_sums[row]], b->row_volts[row]);
}

/* Checks the arc from a variable to a col with a certain value.
   Returns true if the value doesn't break arc consistency. */
bool check_col_arc(struct board *b, int var, int val)
{
	int i, col = COL(var);

	/* GET_BIT(combos[i], j) == 1 iff there is a
	   valid set of assignments to the vars in the
	   col that sums up to i and has j voltorbs.

	   Currently the array is empty, so we
	   construct it in the for loop below. */
	uint32_t combos[16];

	memset(combos + 1, 0, 60);
	combos[0] = 1U;

	if (b->assigns[var] != -1) {
		print_assigns(b);
		printf("Error: checking arc consistency on already-assigned variable %d.\n", var);
		exit(0);
	}

	/* Temporarily assign the variable to the tentative value... */
	b->assigns[var] = val;

	/* Get the possible sums/their voltorbs... */
	for (i = col; i < 25; i += 5)
		incorp_dom(b, i, combos);

	/* And reset the variable to not being assigned. */
	b->assigns[var] = -1;

	/* If there's a combination that has the correct col_sum and col_volts, it's valid. */
	return GET_BIT(combos[b->col_sums[col]], b->col_volts[col]);
}

/* Incorporates the possible values for var
   into a rolling sum stored in array c. */
void incorp_dom(struct board *b, int var, uint32_t *c)
{
	int i, asgmt;
	uint32_t dom, temp[16];

	/* If the value has already been
	   assigned, it only has one value,
	   so we can skip the work below. */
	if (b->assigns[var] >= 0) {

		asgmt = b->assigns[var];

		/* If the assignment is a zero, just left-shift
		   the possible numbers of voltorbs by 1. */
		if (asgmt == 0) {
			for (i = 0; i < 16; i++)
				c[i] <<= 1;
			return;
		}
		
		/* Otherwise, it's the same combos, just
		   shifted up asgmt places in the array. */
		for (i = 15; i >= asgmt; i--)
			c[i] = c[i - asgmt];

		/* Filling in 0s for the bottom section, since
		   we can never reach these values. */
		for (; i >= 0; i--)
			c[i] = 0;

		return;
	}

	/* Otherwise, for every value of var, incorporate it as a
	   possibility by bitwise-ORing all the combinations together. */
	dom = b->domains[var];

	/* The only exception is 0, which we just set temp
	   to, skipping unnecessary bitwise-ORs. */
	if (GET_BIT(dom, 0)) {
		for (i = 0; i < 16; i++)
			temp[i] = c[i] << 1;
	} else {
		memset(temp, 0, 64);
	}

	if (GET_BIT(dom, 1)) {
		for (i = 15; i >= 1; i--)
			temp[i] |= c[i - 1];
	}

	if (GET_BIT(dom, 2)) {
		for (i = 15; i >= 2; i--)
			temp[i] |= c[i - 2];
	}

	if (GET_BIT(dom, 3)) {
		for (i = 15; i >= 3; i--)
			temp[i] |= c[i - 3];
	}

	/* Finally, copy the new combinations into
	   the caller's buffer. */
	memcpy(c, temp, 64);
}

/* Gets the partition constraints for the xth row of the
   given board. Assumes it's passed a non-full row. */
uint32_t row_part_cons(struct board *b, int x)
{
	int i, asgmt;
	uint32_t rv;
	int rem_sum, rem_volts, rem_vars;

	/* Compute data about the unassigned vars in the row.
	   rem_sum => sum of unassigned vars.
	   rem_volts => number of 0s (voltorbs) in unassigned vars.
	   rem_vars => number of unassigned vars, not including voltorbs. */
	rem_sum = b->row_sums[x];
	rem_volts = b->row_volts[x];
	rem_vars = 5 - rem_volts;

	/* Check the assignment status of the vars in the row. */
	for (i = 5 * x; i < 5 * x + 5; i++) {
		asgmt = b->assigns[i];

		/* If it's a voltorb, decrement the number
		   of remaining voltorbs. */
		if (asgmt == 0) {
			rem_volts--;

		/* If it's assigned and nonzero, subtract
		   the value from the remaining total. */
		} else if (asgmt > 0) {
			rem_vars--;
			rem_sum -= asgmt;
		}
	}

	if (rem_vars == 0 && rem_sum != 0) {
		printf("Invalid constraints for row %d. Terminating.\n", x);
		exit(0);
	}

	/* Look up what values the remaining vars can take. */
	rv = GET_PART_CONS(rem_vars, rem_sum);

	/* If there are any voltorbs remaining, any of the
	   remaining vars could still be voltorbs. */
	return (rem_volts > 0 ? SET_BIT(rv, 0) : rv);
}

/* Gets the partition constraints for the xth col of the given board.
   Assumes it's passed a non-full column. */
uint32_t col_part_cons(struct board *b, int x)
{
	int i, asgmt;
	uint32_t rv;
	int rem_sum, rem_volts, rem_vars;

	/* Compute data about the unassigned vars in the col.
	   rem_sum => sum of unassigned vars.
	   rem_volts => number of 0s (voltorbs) in unassigned vars.
	   rem_vars => number of unassigned vars, not including voltorbs. */
	rem_sum = b->col_sums[x];
	rem_volts = b->col_volts[x];
	rem_vars = 5 - rem_volts;

	/* Check the assignment status of the vars in the col. */
	for (i = x; i < 25; i += 5) {
		asgmt = b->assigns[i];

		/* If it's a voltorb, decrement the number
		   of remaining voltorbs. */
		if (asgmt == 0) {
			rem_volts--;

		/* If it's assigned and nonzero, subtract
		   the value from the remaining total. */
		} else if (asgmt > 0) {
			rem_vars--;
			rem_sum -= asgmt;
		}
	}

	if (rem_vars == 0 && rem_sum != 0) {
		printf("Invalid constraints for col %d. Terminating.\n", x);
		exit(0);
	}

	/* Look up what values the remaining vars can take. */
	rv = GET_PART_CONS(rem_vars, rem_sum);

	/* If there are any voltorbs remaining, any of the
	   remaining vars could still be voltorbs. */
	return (rem_volts > 0 ? SET_BIT(rv, 0) : rv);
}

/* Changes a variable's domain to dom, and
   marks any constraints that must be re-checked.
   Returns true if the change resulted in an assignment. */
bool change_var_dom(struct board *b, int var, uint32_t dom)
{
	int i;

	/* If the new domain is the same, no need to do anything. */
	if (dom == b->domains[var])
		return false;

	/* Otherwise, set the variable's domain to the new domain.
	   This also means we need to re-check arc consistency for
	   all non-assigned variables in the row/column. */
	b->domains[var] = dom;
	for (i = 0; i < 5; i++) {
		if (b->assigns[VAR(ROW(var), i)] == -1)
			b->row_arc[VAR(ROW(var), i)] = false;

		if (b->assigns[VAR(i, COL(var))] == -1)
			b->col_arc[VAR(i, COL(var))] = false;
	}

	/* Finally, if it's one-hot, assign the
	   variable and invalidate its row/column. */
	if (dom == 1U) {
		assign_var_doms(b, var, 0);
		return true;
	} else if (dom == 2U) {
		assign_var_doms(b, var, 1);
		return true;
	} else if (dom == 4U) {
		assign_var_doms(b, var, 2);
		return true;
	} else if (dom == 8U) {
		assign_var_doms(b, var, 3);
		return true;
	}

	return false;
}

/* Assigns a value to a variable and invalidates the row/col
   part constraints, unless the row has just been filled.
   Also ensures we no longer check it for arc consistency. */
void assign_var_doms(struct board *b, int var, int val)
{
	b->assigns[var] = val;

	b->row_cons[ROW(var)] = (++(b->row_asgmt[ROW(var)]) != 5) ? false : true;
	b->col_cons[COL(var)] = (++(b->col_asgmt[COL(var)]) != 5) ? false : true;

	b->row_arc[var] = true;
	b->col_arc[var] = true;
}

/* Assigns a value to a variable and sets mask[i]
   to true for every i for which the ith
   solution disagrees with the assignment. */
void assign_var_mask(struct board *b, int var, int val, bool *mask)
{
	int i, len = b->solns_len, *solns = b->solns;

	/* Remove all solutions that disagree with the users input. */
	for (i = 0; i < len; i++) {

		/* Skip all already invalid entries. */
		if (mask[i])
			continue;

		/* If the ith solution disagrees with the
		   assignment, invalidate it in the mask. */
		if (SOLN_ASGMT(solns, i, var) != val)
			mask[i] = true;
	}

}

/* Adds a new solution to the linked list whose tail
   is in **last. Returns true if allocation fails. */
bool ll_append(int *s, struct soln_elem **last)
{
	struct soln_elem *new_soln = malloc(sizeof(struct soln_elem));
	if (!new_soln)
		return true;

	/* Set up the new tail of the list. */
	new_soln->next = NULL;
	memcpy(new_soln->assigns, s, 25 * sizeof(int));

	/* Update the prior last element and point last at the new tail. */
	(*last)->next = new_soln;
	*last = new_soln;

	return false;
}

/* Prints out the constraints on the board (row/column sums/voltorbs) */
void print_cons(struct board *b)
{
	int i, *rs = b->row_sums, *rv = b->row_volts, *cs = b->col_sums, *cv = b->col_volts;
	printf("\nConstraints:\n");

	printf("\nRows:\n");
	for (i = 0; i < 5; i++)
		printf("%02d ", rs[i]);
	printf("\n");
	for (i = 0; i < 5; i++)
		printf(" %d ", rv[i]);

	printf("\n\nCols:\n");
	for (i = 0; i < 5; i++)
		printf("%02d ", cs[i]);
	printf("\n");
	for (i = 0; i < 5; i++)
		printf(" %d ", cv[i]);
}

/* Prints out a grid representation of the assignments of a board. */
void print_assigns(struct board *b)
{
	int i, j, *asgmt = b->assigns;

	printf("\nAssignments:\n\n");
	for (i = 0; i < 5; i++) {
		for (j = 0; j < 5; j++) {
			
			/* Print out X's for unassigned vars, and
			   ?'s for invalidly assigned vars. */
			if (asgmt[VAR(i, j)] == -1)
				printf("X ");
			else if (asgmt[VAR(i, j)] <= 3)
				printf("%d ", asgmt[VAR(i, j)]);
			else
				printf("? ");
		}

		printf("\n\n");
	}

	printf("\n\n");
}

/* Prints out the domains of all variables for the given board.
   These may not be accurate past the arc-consistency phase. */
void print_doms(struct board *b)
{
	int i, j;
	uint32_t *doms = b->domains;

	printf("\nDomains:\n\n");
	for (i = 0; i < 5; i++) {
		for (j = 0; j < 5; j++)
			printf("%u ", doms[VAR(i, j)]);

		printf("\n\n");
	}

	printf("\n\n");
}

/* Prints out the possible solutions to a
   board, if there are any populated. */
void print_solns(struct board *b)
{
	int i, j, k, l = 0, len = b->solns_len, rem = 0, *solns = b->solns;
	bool *mask = b->mask;

	printf("Size of underlying array: 25 * %d = %d\n\n", len, 25 * len);

	for (i = 0; i < len; i++) {

		/* If it's no longer valid, skip it. */
		if (mask[i])
			continue;

		rem++;

		/* Otherwise, print out a grid for the board. */
		printf("Solution %d:\n\n", ++l);
		for (j = 0; j < 5; j++) {
			for (k = 0; k < 5; k++)
				printf("%d ", SOLN_ASGMT(solns, i, VAR(j, k)));

			printf("\n\n");
		}

		printf("\n\n");
	}

	printf("Remaining solutions: %d\n\n", rem);
}

/* Prints out the probability distribution of all
   variables that aren't already assigned. */
void print_probs(double *probs)
{
	int i, j;

	for (i = 0; i < 25; i++) {

		/* If it's already assigned, skip it. */
		if (PROB(probs, i, 0) == 1.0 ||
				PROB(probs, i, 1) == 1.0 ||
				PROB(probs, i, 2) == 1.0 ||
				PROB(probs, i, 3) == 1.0)
			continue;

		/* Otherwise, print out the probabilities for each of its values. */
		printf("Probabilities for variable %d:\n", i);
		for (j = 0; j < 4; j++)
			printf("%d: %01.02f     ", j, PROB(probs, i, j));

		printf("\n\n");
	}
}

