/* TODO */
/* Incorporate the value of the next game. This will
   require some info about expected winrates, etc. */
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>

/* Set this flag to see (large) log output. */
#define LOG 0

/* Maximum depth of tree search. Larger makes better decisions, but
   runtime is exponential wrt. this argument. */
#define MAX_DEPTH 4

/* Bitmap corresponding to a full domain for a variable. */
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

/* Macro for the evaluation function so we can easily change it. */
#define EVAL(b, m) score(b, m)

/* Error codes. */
enum errs {
	/* Debug codes. */
	NOT_YET_IMP = 1,
	/* Runtime error codes. */
	BAD_ARGS,
	OOM
};

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
   the first four arrays. Early on in the exploration process, our
   knowledge about assignments is encoded in domains. After
   enumeration, this array becomes deprecated, and instead our
   knowledge is stored in solns, a flattened array that contains 
   all possible solutions to the constraints given what we know.
   As we gain more info, rather than re-alloc the array, we use
   mask to indicate which solutions are no longer valid. */
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

	/* A bitmap corresponding to the remaining domains of
	   each variable. Specifically, the nth bit is set
	   if n is a valid value for the variable.
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
	   we just set mask[i] to invalidate it. */
	int *solns, solns_len;
	bool *mask;

	/* An array holding the probabilities of each
	   value occuring in a certain variable. */
	double probs[100];
};

/* Algorithmic functions. */
static void init_board(struct board *b, int *rs, int *rv, int *cs, int *cv);
static int solve_board(struct board *b);
static int enf_part(struct board *b, bool *rd);
static bool enf_arcs(struct board *b);
static bool pfree_doms(struct board *b);
static int enum_boards(struct board *b);
static void populate_probs(struct board *b, bool *mask, double *probs);
static bool pfree_mask(struct board *b);
static int pbest_act(struct board *b, bool *rd);
static int row_part_cons(struct board *b, int x, uint32_t *rd);
static int col_part_cons(struct board *b, int x, uint32_t *rd);
static bool change_var_dom(struct board *b, int var, uint32_t dom);
static bool check_row_arc(struct board *b, int var, int val);
static bool check_col_arc(struct board *b, int var, int val);
static int r_enum_boards(struct board *b, struct soln_elem **last, int depth);
static void assign_var_mask(struct board *b, int var, int val, bool *mask);
static void update_asgmt(struct board *b);
static bool is_solved(struct board *b, bool *mask);
static int tree_search(struct board *b, bool *m_in, int var, int val, int depth, double *rd);
static void assign_var_doms(struct board *b, int var, int val);
static void incorp_dom(struct board *b, int var, uint32_t *c);
static int ll_append(int *s, struct soln_elem **last);
static int score(struct board *b, bool *mask);
static void print_assigns(struct board *b);
static void print_prompt(struct board *b, int x);
static void print_error(enum errs code);

#if LOG
/* Log functions. */
static void print_doms(struct board *b);
static void print_solns(struct board *b);
static void print_probs(double *probs);
#endif

/* Main. */
int main()
{
	int errno;
	struct board _b;
	struct board *b = &_b;
#if 0
	/* 3-rooks. Close to maximum uncertainty
	   where quitting is still unadvisable. */
	int rs[] = {2, 2, 2, 0, 1};
	int rv[] = {4, 4, 4, 5, 4};
	int cs[] = {2, 2, 2, 0, 1};
	int cv[] = {4, 4, 4, 5, 4};

	/* 4-rooks. Quitting is advisable
	   after clicking the known 1 tile. */
	int rs[] = {2, 2, 2, 2, 1};
	int rv[] = {4, 4, 4, 4, 4};
	int cs[] = {2, 2, 2, 2, 1};
	int cv[] = {4, 4, 4, 4, 4};

	/* Same as above, but since our
	   beginning score is 0 it's worth
	   trying a tile. */
	int rs[] = {2, 2, 2, 2, 0};
	int rv[] = {4, 4, 4, 4, 5};
	int cs[] = {2, 2, 2, 2, 0};
	int cv[] = {4, 4, 4, 4, 5};
#endif
	int rs[] = {4, 5, 3, 7, 7};
	int rv[] = {2, 2, 4, 1, 1};
	int cs[] = {2, 7, 5, 8, 4};
	int cv[] = {3, 1, 2, 1, 3};

	init_board(b, rs, rv, cs, cv);

	errno = solve_board(b);
	if (errno) {
		print_error(errno);
		printf("\nBoard state at failure:\n");
		print_assigns(b);
	}

	free(b->mask);
	free(b->solns);
	return 0;
}

/* Initializes a board with the given row and column constraints. */
static void init_board(struct board *b, int *rs, int *rv, int *cs, int *cv)
{
	int i, *asgmt = b->assigns;
	uint32_t *dom = b->domains;

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
static int solve_board(struct board *b)
{
	int errno;
	bool cont;

	/* Reduce using arc consistency, since enum_boards
	   expects an already arc-reduced board. */
	while (1) {

		/* First, reduce as much as possible using partition constraints. */
		errno = enf_part(b, &cont);
		if (errno) {
			printf("ERROR: Failed to enforce partition constraints: %d.\n", __LINE__);
			return errno;
		}

		if (cont)
			continue;

		/* Reduce as much as possible using arc consistency. Finally, prompt
		   the user for any safe tiles. When none of the three steps in this
		   loop produce a new assignment, break out and continue. */
		if (!enf_arcs(b) && !pfree_doms(b))
			break;
	}
#if LOG
	printf("\nDomains after one round of arc consistency:\n");
	print_doms(b);
#endif
	/* Create a list of all possible boards given the constraints. */
	errno = enum_boards(b);
	if (errno) {
		printf("ERROR: Failed to enumerate all possible solutions to the board: %d.\n", __LINE__);
		return errno;
	}

	/* Populate the probabilities so that we can pass
	   this board into pfree_mask(). */
	populate_probs(b, b->mask, b->probs);
#if LOG
	print_solns(b);
	print_assigns(b);
	printf("Probabilities of unassigned tiles:\n");
	print_probs(b->probs);
#endif
	/* Prompt the user to reveal all tiles that are safe,
	   until there are no more safe tiles. Then, do the
	   best action as given by our game theory engine.
	   Repeat until the game is over. */
	while (1) {
		if (pfree_mask(b))
			continue;

		/* Prompt the user to reveal a tile. */
		errno = pbest_act(b, &cont);
		if (errno) {
			free(b->solns);
			free(b->mask);
			printf("ERROR: Failed to prompt the user for the next best action: %d.\n", __LINE__);
			return errno;
		}

		/* If the game is over, we're done. */
		if (!cont)
			return 0;
	}
}

/* Changes b->domains to enforce partition constraints.
   Sets *rd to true if an assignment was made. */
static int enf_part(struct board *b, bool *rd)
{
	int errno, i, j;
	uint32_t cons;
	bool rv = false;

	/* Iterates through the rows of the board. */
	for (i = 0; i < 5; i++) {

		/* If the row has already been enforced, skip it. */
		if (b->row_cons[i])
			continue;

		b->row_cons[i] = true;

		/* Look up the constraints on the variables. */
		errno = row_part_cons(b, i, &cons);
		if (errno) {
			printf("ERROR: Could not obtain partition constraints for row %d: %d.\n", i, __LINE__);
			return errno;
		}

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
		errno = col_part_cons(b, i, &cons);
		if (errno) {
			printf("ERROR: Could not obtain partition constraints for col %d: %d.\n", i, __LINE__);
			return errno;
		}

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

	*rd = rv;
	return 0;
}

/* Enforces arc consistency on all variables in board b.
   Returns as soon as it makes an assignment so enf_part
   can take over and enforce partitioning before continuing. */
static bool enf_arcs(struct board *b)
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
			if ((!b->row_arc[i] && !check_row_arc(b, i, j)) ||
					(!b->col_arc[i] && !check_col_arc(b, i, j)))
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

/* Prompts the user for the true value of a safe tile on the board.
   This version operates on the domains framework, to reduce the
   amount of memory allocated for the solutions list. */
static bool pfree_doms(struct board *b)
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
	print_prompt(b, free_var);
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
static int enum_boards(struct board *b)
{
	int errno, i, num, *solns;
	struct soln_elem ll_sent, *soln_iter, *soln_next, *last;
	bool *mask;

	/* Initialize the linked list sentinel node. */
	ll_sent.next = NULL;
	last = &(ll_sent);

	/* Populates the linked list with all solutions to b. */
	errno = r_enum_boards(b, &last, 1);
	num = b->solns_len;
	if (errno) {
		/* Free whatever portions of the LL we allocated. */
		for (i = 0, soln_iter = ll_sent.next; i < num; i++, soln_iter = soln_next) {
			soln_next = soln_iter->next;
			free(soln_iter);
		}

		printf("ERROR: Failed to enumerate boards with backtracking DFS: %d.\n", __LINE__);
		return errno;
	}

	/* Malloc the (usually large) set of solutions and their mask. */
	solns = malloc(25 * num * sizeof(int));
	if (!solns) {
		/* Free the linked list. */
		for (i = 0, soln_iter = ll_sent.next; i < num; i++, soln_iter = soln_next) {
			soln_next = soln_iter->next;
			free(soln_iter);
		}

		printf("ERROR: Failed to allocate solution buffer: %d.\n", __LINE__);
		return OOM;
	}

	mask = malloc(num * sizeof(bool));
	if (!mask) {
		/* Free the linked list. */
		for (i = 0, soln_iter = ll_sent.next; i < num; i++, soln_iter = soln_next) {
			soln_next = soln_iter->next;
			free(soln_iter);
		}

		free(solns);
		printf("ERROR: Failed to allocate mask buffer: %d.\n", __LINE__);
		return OOM;
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

	b->mask = mask;
	b->solns = solns;
	return 0;
}

/* Populates an array with the probability of each value
   of every variable given a mask for the board b. */
static void populate_probs(struct board *b, bool *mask, double *probs)
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

/* If there's a safe variable on the board, prompt the user
   for its value. Operates on the solutions enum/mask framework. */
static bool pfree_mask(struct board *b)
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
	printf("\n##################################\n");
	print_prompt(b, free_var);
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

/* Finds the best move on the board by using MAX_DEPTH game tree search, and
   prompts the user to take this action. Sets *rd to true if the game continues,
   and false if it ends (due to quitting, winning, or clicking a voltorb). */
static int pbest_act(struct board *b, bool *rd)
{
	int errno, i, input, best_a, *asgmt = b->assigns;
	double best = -1.0, best_pzero = 1.1, cand, p_one, p_two, p_three, sub_val, *probs = b->probs;
	uint32_t dom = 0;
#if LOG
	int j;
	double q_vals[25];
#endif
	/* If the board is already solved, we're done! */
	if (is_solved(b, b->mask)) {
		print_assigns(b);
		printf("All 2s and 3s have been found - the board is solved!\n\n");
		*rd = false;
		return 0;
	}

	/* Otherwise, find the next action. First, compute the
	   expected values of revealing all the unknown tiles. */
	for (i = 0; i < 25; i++) {
#if LOG
		q_vals[i] = 0.0;
#endif
		/* If the variable is already assigned, we
		   don't gain anything by revealing it. */
		if (asgmt[i] != -1)
			continue;
#if LOG
		printf("Assigning var %d.\n\n", i);
#endif
		/* Otherwise, calculate its expected value. */
		cand = 0;
		p_one = PROB(probs, i, 1);
		p_two = PROB(probs, i, 2);
		p_three = PROB(probs, i, 3);

		/* If the probability of a certain value is > 0.0, add the value of
		   the subtree multiplied by its probability to the expected value
		   of the action. Note that 0 always results in 0, so we can skip it. */
		if (p_one > 0.0) {
			errno = tree_search(b, b->mask, i, 1, 1, &sub_val);
			if (errno) {
				printf("ERROR: Failed to compute value of subtree 1 -> %d at the root: %d.\n", i, __LINE__);
				return errno;
			}

			cand += p_one * sub_val;
		}

		/* Same as above, but for 2. */
		if (p_two > 0.0) {
			errno = tree_search(b, b->mask, i, 2, 1, &sub_val);
			if (errno) {
				printf("ERROR: Failed to compute value of subtree 2 -> %d at the root: %d.\n", i, __LINE__);
				return errno;
			}

			cand += p_two * sub_val;
		}

		/* Again, but for 3. */
		if (p_three > 0.0) {
			errno = tree_search(b, b->mask, i, 3, 1, &sub_val);
			if (errno) {
				printf("ERROR: Failed to compute value of subtree 3 -> %d at the root: %d.\n", i, __LINE__);
				return errno;
			}

			cand += p_three * sub_val;
		}
#if LOG
		q_vals[i] = cand;
		printf("Value of clicking on %d at the root is %01.02f.\n\n", i, cand);
#endif
		/* If it's better, set this as the best action so far. */
		if (cand > best || (cand == best && PROB(probs, i, 0) < best_pzero)) {
			best = cand;
			best_pzero = PROB(probs, i, 0);
			best_a = i;
		}
	}
#if LOG
	printf("Q-Values:\n\n");
	for (i = 0; i < 5; i++) {
		for (j = 0; j < 5; j++)
			printf("%01.02f ", q_vals[VAR(i, j)]);
		printf("\n");
	}
#endif
	/* If quitting yields a higher score than revealing any
	   particular tile, prompt the user to quit and end the game. */
	if (score(b, b->mask) > best) {
		print_assigns(b);
		printf("It is not worthwhile to continue. Maximize your score with the above known tiles and quit the game.\n\n");
		*rd = false;
		return 0;
	}

	printf("\n##################################\n");
	print_prompt(b, best_a);
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

	printf("\nIt has a %01.02f%% chance to be a Voltorb.", 100 * PROB(probs, best_a, 0));

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

	/* Updates the mask of the board, removing all
	   solns that disagree with input. */
	assign_var_mask(b, best_a, input, b->mask);

	/* Officially assign the value to the variable. */
	asgmt[best_a] = input;

	/* Recalculate the probabilities given the new list of boards. */
	populate_probs(b, b->mask, b->probs);
	update_asgmt(b);

	/* If they clicked on a voltorb, the game is over. */
	if (input == 0) {
		printf("\n##################################\n");
		print_assigns(b);
		printf("Sorry, you got unlucky. =( Game Over.\n\n");
		*rd = false;
		return 0;
	}
#if LOG
	print_assigns(b);
	printf("Probabilities of unassigned tiles:\n");
	print_probs(b->probs);
#endif
	*rd = true;
	return 0;
}

/* Gets the partition constraints for the xth row of the
   given board. Assumes it's passed a non-full row. */
static int row_part_cons(struct board *b, int x, uint32_t *rd)
{
	int i, asgmt, rem_sum, rem_volts, rem_vars;
	uint32_t rv;

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

	/* Sanity test for invalid constraints. */
	if (rem_vars == 0 && rem_sum != 0) {
		printf("ERROR: Invalid constraints for row %d: %d.\n", x, __LINE__);
		return BAD_ARGS;
	}

	/* Look up what values the remaining vars can take. */
	rv = GET_PART_CONS(rem_vars, rem_sum);

	/* If there are any voltorbs remaining, any of the
	   remaining vars could still be voltorbs. */
	if (rem_volts > 0)
		rv = SET_BIT(rv, 0);

	*rd = rv;
	return 0;
}

/* Gets the partition constraints for the xth col of the given board.
   Assumes it's passed a non-full column. */
static int col_part_cons(struct board *b, int x, uint32_t *rd)
{
	int i, asgmt, rem_sum, rem_volts, rem_vars;
	uint32_t rv;

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

	/* Sanity test for invalid constraints. */
	if (rem_vars == 0 && rem_sum != 0) {
		printf("ERROR: Invalid constraints for col %d: %d.\n", x, __LINE__);
		return BAD_ARGS;
	}

	/* Look up what values the remaining vars can take. */
	rv = GET_PART_CONS(rem_vars, rem_sum);

	/* If there are any voltorbs remaining, any of the
	   remaining vars could still be voltorbs. */
	if (rem_volts > 0)
		rv = SET_BIT(rv, 0);

	*rd = rv;
	return 0;
}

/* Changes a variable's domain to dom, and marks
   any constraints that must be re-enforced.
   Returns true if the change resulted in an assignment. */
static bool change_var_dom(struct board *b, int var, uint32_t dom)
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

/* Checks the arc from a variable to a row with a certain value.
   Returns true if the value doesn't break arc consistency. */
static bool check_row_arc(struct board *b, int var, int val)
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
static bool check_col_arc(struct board *b, int var, int val)
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

/* Populates the given linked list with all solutions that satisfy
   the board's constraints using backtracking search with arc
   consistency. Assumes it's passed an already arc-reduced board. */
static int r_enum_boards(struct board *b, struct soln_elem **last, int depth)
{
	int errno, i, j, asgmts[25], card, min_card = 5, min_asgmt = -1;
	uint32_t doms[25], dom;
	bool cont, valid;

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
		errno = ll_append(asgmts, last);

		/* If an error occurs, report the depth that we failed to ll_append. */
		if (errno)
			printf("ERROR: Failed to append a solution to the LL of solutions at depth %d: %d.\n", depth, __LINE__);

		return errno;
	}

	/* Otherwise, try on each value for the most-constrained
	   variable, and recursively call r_enum_boards. */
	dom = doms[min_asgmt];
	for (i = 0; i < 4; i++) {

		/* Skip invalid assignments. */
		if (!GET_BIT(dom, i))
			continue;

		/* Assign the variable to i by restricting its domain.
		   This also invalidates arc-consistency. */
		change_var_dom(b, min_asgmt, 1U << i);

		/* Re-enforce arc consistency... */
		while (1) {

			/* First, reduce as much as possible using partition constraints. */
			errno = enf_part(b, &cont);
			if (errno) {
				printf("ERROR: Failed to enforce partition constraints: %d.\n", __LINE__);
				return errno;
			}

			if (cont)
				continue;

			/* Reduce as much as possible using arc consistency.
			   Unlike in solve_board, it makes no sense to prompt the user here. */
			if (!enf_arcs(b))
				break;
		}

		/* If any variables have empty domains, this is an invalid assignment. */
		valid = true;
		for (j = 0; j < 25; j++) {
			if (b->domains[j] == 0)
				valid = false;
		}

		/* Recursively call ourselves, but only if it's a valid assignment. */
		if (valid) {
			errno = r_enum_boards(b, last, depth + 1);

			/* If an error occurs, pass it up. */
			if (errno) {
				printf("ERROR: Failed to enumerate subtree after %d -> %d at depth %d: %d.\n", i, min_asgmt, depth, __LINE__);
				return errno;
			}
		}

		/* And restore the previous state
		   of b for backtracking. */
		memcpy(b->assigns, asgmts, 25 * sizeof(int));
		memcpy(b->domains, doms, 100);
	}

	return 0;
}

/* Assigns a value to a variable and sets mask[i] to true for every
   i for which the ith solution disagrees with the assignment. */
static void assign_var_mask(struct board *b, int var, int val, bool *mask)
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

/* Returns whether a given board has been solved - that is,
   whether or not all the hidden 2s and 3s have been found. */
static bool is_solved(struct board *b, bool *mask)
{
	int i, j, len = b->solns_len, *solns = b->solns,
	    asgmt[25], soln_asgmt;

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

/* Say we have a board defined by b and m_in. This function either
   returns the expected value of the best action, or, if the board is
   solved (or we've hit MAX_DEPTH), the evaluation of the board. */
static int tree_search(struct board *b, bool *m_in, int var, int val, int depth, double *rd)
{
	int errno, i, len = b->solns_len;
	double rv = -1.0, cand, p_one, p_two, p_three, sub_val, probs[100];
	bool *mask = malloc(len * sizeof(bool));
#if LOG
	int _;
#endif
	if (!mask) {
		printf("ERROR: Ran out of memory trying to malloc a mask for %d -> %d at depth %d: %d\n", val, var, depth, __LINE__);
		return OOM;
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

		if (depth == MAX_DEPTH)
			printf("Reached maximum recursion depth! Value of board is %01.02f.\n\n", (double) EVAL(b, mask));
		else
			printf("Found a solution! Value of board is %01.02f.\n\n", (double) EVAL(b, mask));
#endif
		*rd = EVAL(b, mask);
		free(mask);
		return 0;
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
		printf("Assigning var %d.\n\n", i);
#endif
		/* Otherwise, calculate its expected value. */
		cand = 0;
		p_one = PROB(probs, i, 1);
		p_two = PROB(probs, i, 2);
		p_three = PROB(probs, i, 3);

		/* If the probability of a certain value is > 0.0, add the value of
		   the subtree multiplied by its probability to the expected value
		   of the action. Note that 0 always results in 0, so we can skip it. */
		if (p_one > 0.0) {
			errno = tree_search(b, mask, i, 1, depth + 1, &sub_val);
			if (errno) {
				printf("ERROR: Failed to compute value of subtree 1 -> %d at depth %d: %d.\n", i, depth, __LINE__);
				free(mask);
				return errno;
			}

			cand += p_one * sub_val;
		}

		/* Same as above, but for 2. */
		if (p_two > 0.0) {
			errno = tree_search(b, mask, i, 2, depth + 1, &sub_val);
			if (errno) {
				printf("ERROR: Failed to compute value of subtree 2 -> %d at depth %d: %d.\n", i, depth, __LINE__);
				free(mask);
				return errno;
			}

			cand += p_two * sub_val;
		}

		/* Again, but for 3. */
		if (p_three > 0.0) {
			errno = tree_search(b, mask, i, 3, depth + 1, &sub_val);
			if (errno) {
				printf("ERROR: Failed to compute value of subtree 3 -> %d at depth %d: %d.\n", i, depth, __LINE__);
				free(mask);
				return errno;
			}

			cand += p_three * sub_val;
		}
#if LOG
		for (_ = 0; _ < depth; _++)
			printf("\t");
		printf("Value of clicking on %d at depth %d is %01.02f.\n\n", i, depth, cand);
#endif
		/* If the expected value for this action is the best we've found so far, remember it. */
		if (cand > rv)
			rv = cand;
	}

	/* Check if quitting now would be the better option. */
	if (score(b, mask) > rv)
		rv = score(b, mask);

	free(mask);
	*rd = rv;
	return 0;
}

/* Assigns a value to a variable and invalidates the row/col
   part constraints, unless the row has just been filled.
   Also ensures we no longer check it for arc consistency. */
static void assign_var_doms(struct board *b, int var, int val)
{
	b->assigns[var] = val;

	b->row_cons[ROW(var)] = (++(b->row_asgmt[ROW(var)]) != 5) ? false : true;
	b->col_cons[COL(var)] = (++(b->col_asgmt[COL(var)]) != 5) ? false : true;

	b->row_arc[var] = true;
	b->col_arc[var] = true;
}

/* Incorporates the possible values for var
   into a rolling sum stored in array c. */
static void incorp_dom(struct board *b, int var, uint32_t *c)
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

/* Adds a new solution to the linked list whose tail
   is in **last. Spawns an error if allocation fails. */
static int ll_append(int *s, struct soln_elem **last)
{
	struct soln_elem *new_soln = malloc(sizeof(struct soln_elem));
	if (!new_soln) {
		printf("ERROR: Failed to allocate linked-list element: %d.\n", __LINE__);
		return OOM;
	}

	/* Set up the new tail of the list. */
	new_soln->next = NULL;
	memcpy(new_soln->assigns, s, 25 * sizeof(int));

	/* Update the prior last element and point last at the new tail. */
	(*last)->next = new_soln;
	*last = new_soln;

	return 0;
}

/* Returns the maximum score we can safely obtain on b knowing mask.
   This is the score one would receive if they quit at this point. */
static int score(struct board *b, bool *mask)
{
	int i, j, rv = 1, len = b->solns_len, *solns = b->solns;
	uint32_t doms[25];
	bool zero_score = true;

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
		if (doms[i] == 2U) {
			zero_score = false;
		} else if (doms[i] == 4U) {
			zero_score = false;
			rv *= 2;
		} else if (doms[i] == 8U) {
			zero_score = false;
			rv *= 3;
		}
	}

	return zero_score ? 0 : rv;
}

/* Prints out a grid representation of the
   known variable assignments for  a board. */
static void print_assigns(struct board *b)
{
	int i, j, *asgmt = b->assigns, *rs = b->row_sums, *rv = b->row_volts,
	    *cs = b->col_sums, *cv = b->col_volts;

	printf("\nAssignments:\n\n");
	for (i = 0; i < 5; i++) {

		printf("\n");
		for (j = 0; j < 5; j++) {
			
			/* Print out "X"s for unassigned vars. */
			if (asgmt[VAR(i, j)] == -1)
				printf(" X ");
			else
				printf(" %d ", asgmt[VAR(i, j)]);
		}

		printf("   %02d\n                  %d\n\n", rs[i], rv[i]);
	}

	printf("%02d %02d %02d %02d %02d\n %d  %d  %d  %d  %d\n\n",
			cs[0], cs[1], cs[2], cs[3], cs[4],
			cv[0], cv[1], cv[2], cv[3], cv[4]);
}

/* Does the same thing as print_assigns, but with a convenient
   ? instead of an X for the given tile. */
static void print_prompt(struct board *b, int x)
{
	int i, j, *asgmt = b->assigns, *rs = b->row_sums, *rv = b->row_volts,
	    *cs = b->col_sums, *cv = b->col_volts;

	printf("\nAssignments:\n\n");
	for (i = 0; i < 5; i++) {

		printf("\n");
		for (j = 0; j < 5; j++) {

			if (VAR(i, j) == x) {
				printf(" ? ");
				continue;
			}
			
			/* Print out "X"s for unassigned vars. */
			if (asgmt[VAR(i, j)] == -1)
				printf(" X ");
			else
				printf(" %d ", asgmt[VAR(i, j)]);
		}

		printf("   %02d\n                  %d\n\n", rs[i], rv[i]);
	}

	printf("%02d %02d %02d %02d %02d\n %d  %d  %d  %d  %d\n\n",
			cs[0], cs[1], cs[2], cs[3], cs[4],
			cv[0], cv[1], cv[2], cv[3], cv[4]);
}

/* Prints out a helpful message corresponding to an error code. */
static void print_error(enum errs code)
{
	switch (code) {
	/* Debug codes. */
	case NOT_YET_IMP:
		printf("Not yet implemented.\n");
		break;
	/* Runtime error codes. */
	case BAD_ARGS:
		printf("Bad arguments.\n");
		break;
	case OOM:
		printf("Ran out of memory.\n");
		break;
	}
}

#if LOG
/* Prints out the domains of all variables for the given board.
   These may not be accurate past the arc-consistency phase. */
static void print_doms(struct board *b)
{
	int i, j, *rs = b->row_sums, *rv = b->row_volts,
	    *cs = b->col_sums, *cv = b->col_volts;

	uint32_t *doms = b->domains;

	printf("\nDomains:\n\n");
	for (i = 0; i < 5; i++) {

		printf("\n");
		for (j = 0; j < 5; j++)
			printf(" %u ", doms[VAR(i, j)]);

		printf("   %02d\n                  %d\n\n", rs[i], rv[i]);
	}

	printf("%02d %02d %02d %02d %02d\n %d  %d  %d  %d  %d\n\n",
			cs[0], cs[1], cs[2], cs[3], cs[4],
			cv[0], cv[1], cv[2], cv[3], cv[4]);
}

/* Prints out the possible solutions to a
   board, if there are any populated. */
static void print_solns(struct board *b)
{
	int i, j, k, l = 0, len = b->solns_len, rem = 0, *solns = b->solns,
	    *rs = b->row_sums, *rv = b->row_volts,
	    *cs = b->col_sums, *cv = b->col_volts;

	bool *mask = b->mask;

	/* If the solns array hasn't been populated, print an error and continue. */
	if (!solns) {
		printf("ERROR: Tried to print out the solutions list of an unpopulated board: %d.\n", __LINE__);
		return;
	}

	for (i = 0; i < len; i++) {

		/* If it's no longer valid, skip it. */
		if (mask[i])
			continue;

		rem++;

		/* Otherwise, print out a grid for the board. */
		printf("##################################\n\nSolution %d:\n\n", ++l);
		for (j = 0; j < 5; j++) {

			printf("\n");
			for (k = 0; k < 5; k++)
				printf(" %u ", SOLN_ASGMT(solns, i, VAR(j, k)));

			printf("   %02d\n                  %d\n\n", rs[j], rv[j]);
		}

		printf("%02d %02d %02d %02d %02d\n %d  %d  %d  %d  %d\n\n",
				cs[0], cs[1], cs[2], cs[3], cs[4],
				cv[0], cv[1], cv[2], cv[3], cv[4]);
	}

	printf("\nSize of underlying array: 25 * %d = %d\nRemaining solutions: %d\n\n", len, 25 * len, rem);
}

/* Prints out the probability distribution of all
   variables that aren't already assigned. */
static void print_probs(double *probs)
{
	int i, j;

	for (i = 0; i < 25; i++) {

		/* If the variable's value is certain, skip it. */
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
#endif

