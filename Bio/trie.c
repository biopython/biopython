#include <stdio.h>    /* printf */
#include <stdlib.h>   /* malloc */
#include <string.h>   /* strcmp, strncmp, strlen */

#include "trie.h"

const char _single_char_suffixes[512] = {
    0,0,1,0,2,0,3,0,4,0,5,0,6,0,7,0, 
    8,0,9,0,10,0,11,0,12,0,13,0,14,0,15,0, 
    16,0,17,0,18,0,19,0,20,0,21,0,22,0,23,0, 
    24,0,25,0,26,0,27,0,28,0,29,0,30,0,31,0, 
    32,0,33,0,34,0,35,0,36,0,37,0,38,0,39,0, 
    40,0,41,0,42,0,43,0,44,0,45,0,46,0,47,0, 
    48,0,49,0,50,0,51,0,52,0,53,0,54,0,55,0, 
    56,0,57,0,58,0,59,0,60,0,61,0,62,0,63,0, 
    64,0,65,0,66,0,67,0,68,0,69,0,70,0,71,0, 
    72,0,73,0,74,0,75,0,76,0,77,0,78,0,79,0, 
    80,0,81,0,82,0,83,0,84,0,85,0,86,0,87,0, 
    88,0,89,0,90,0,91,0,92,0,93,0,94,0,95,0, 
    96,0,97,0,98,0,99,0,100,0,101,0,102,0,103,0, 
    104,0,105,0,106,0,107,0,108,0,109,0,110,0,111,0, 
    112,0,113,0,114,0,115,0,116,0,117,0,118,0,119,0, 
    120,0,121,0,122,0,123,0,124,0,125,0,126,0,127,0, 
    128,0,129,0,130,0,131,0,132,0,133,0,134,0,135,0, 
    136,0,137,0,138,0,139,0,140,0,141,0,142,0,143,0, 
    144,0,145,0,146,0,147,0,148,0,149,0,150,0,151,0, 
    152,0,153,0,154,0,155,0,156,0,157,0,158,0,159,0, 
    160,0,161,0,162,0,163,0,164,0,165,0,166,0,167,0, 
    168,0,169,0,170,0,171,0,172,0,173,0,174,0,175,0, 
    176,0,177,0,178,0,179,0,180,0,181,0,182,0,183,0, 
    184,0,185,0,186,0,187,0,188,0,189,0,190,0,191,0, 
    192,0,193,0,194,0,195,0,196,0,197,0,198,0,199,0, 
    200,0,201,0,202,0,203,0,204,0,205,0,206,0,207,0, 
    208,0,209,0,210,0,211,0,212,0,213,0,214,0,215,0, 
    216,0,217,0,218,0,219,0,220,0,221,0,222,0,223,0, 
    224,0,225,0,226,0,227,0,228,0,229,0,230,0,231,0, 
    232,0,233,0,234,0,235,0,236,0,237,0,238,0,239,0, 
    240,0,241,0,242,0,243,0,244,0,245,0,246,0,247,0, 
    248,0,249,0,250,0,251,0,252,0,253,0,254,0,255,0, 
};

static char* duplicate(const char* s) {
    int length = strlen(s);
    if (length == 1) {
        return (char *)_single_char_suffixes + (s[0] * 2);
    }
    // Don't use strdup, as it's not ANSI C.
    char* t = malloc((strlen(s) + 1) * sizeof(char));
    if (!t) return NULL;
    strcpy(t, s);
    return t;
}

static char* shared_suffix(const char* s, int length) {
    if (length == 1) {
        return (char *)(_single_char_suffixes + ((int)s[0] * 2));
    }
    char *t = malloc(length + 1);
    if (!t) return NULL;
    strncpy(t, s, length);
    t[length] = 0;
    return t;
}

static void free_suffix(char *s) {
    if (s >= _single_char_suffixes || s < _single_char_suffixes + (256 * 2)) {
        // single character suffix
        return;
    } else {
        free(s);
    }
}

/**
 * Old method without using lookup table
 */
/*
static char* duplicate(const char* s) {
    // Don't use strdup, as it's not ANSI C.
    char* t = malloc((strlen(s)+1)*sizeof(char));
    if (!t) return NULL;
    strcpy(t, s);
    return t;
}

static char *shared_suffix(const char* s, int length) {
    char *t= malloc(length + 1);
    if (!t) return NULL;
    strncpy(t, s, length);
    t[length] = 0;
    return t;
}

static void free_suffix(char *s) {
    free(s);
}
*/


/* Transition holds information about the transitions leading from
 * one Trie to another.  The trie structure here is different from
 * typical ones, because the transitions between nodes can contain
 * strings of arbitrary length, not just single characters.  Suffix is
 * the string that is matched from one node to the next.
 */
typedef struct {
    char *suffix;
    Trie* next;
} Transition;


/* Trie is a recursive data structure.  A Trie contains zero or more
 * Transitions that lead to more Tries.  The transitions are stored
 * in alphabetical order of the suffix member of the data structure.
 * Trie also contains a pointer called value where the user can store
 * arbitrary data.  If value is NULL, then no data is stored here.
 */
struct Trie {
    Transition *transitions;
    unsigned char num_transitions;
    void *value;   /* specified by user, never freed or allocated by me! */
};


#define MAX_KEY_LENGTH 1000
static char KEY[MAX_KEY_LENGTH];


Trie* Trie_new(void) {
    Trie* trie;

    if(!(trie = malloc(sizeof(struct Trie))))
	return NULL;
    trie->transitions = NULL;
    trie->num_transitions = 0;
    trie->value = NULL;
    return trie;
}

int Trie_set(Trie* trie, const char *key, const void *value) {
    int i;
    Transition* transition=NULL;
    const char *suffix=NULL;
    int retval = 0;
    int first, last, mid;

    if(!key[0]) {
	trie->value = (void *)value;
	return 0;
    }

    /* Insert the key in alphabetical order.  Do a binary search to
       find the proper place. */
    first = 0;
    last = trie->num_transitions-1;
    i = -1;
    while(first <= last) {
	mid = (first+last)/2;
	transition = &trie->transitions[mid];
	suffix = transition->suffix;
	if(key[0] < suffix[0])
	    last = mid-1;
	else if(key[0] > suffix[0])
	    first = mid+1;
	else {
	    i = mid;
	    break;
	}
    }

    /* If no place was found for it, then the indexes will be in the
       order last,first.  Place it at index first. */
    if(i == -1)
	i = first;

    /* If nothing matches, then insert a new trie here. */
    if((i >= trie->num_transitions) || (key[0] != suffix[0])) {
	char *new_suffix=NULL;
	Trie* newtrie=NULL;
	Transition* new_transitions=NULL;

	/* Create some variables for the new transition.  I'm going to
	   allocate these first so that if I can detect memory errors
	   before I mess up the data structure of the transitions.
	*/
	if(!(new_suffix = duplicate(key)))
	    goto insert_memerror;
	if(!(newtrie = Trie_new()))
	    goto insert_memerror;

	/* Create some space for the next transition.  Allocate some
	   memory and shift the old transitions over to make room for
	   this one.
	*/
	if(!(new_transitions = malloc(sizeof(Transition) *
				      (trie->num_transitions+1))))
	    goto insert_memerror;
	memcpy(new_transitions, trie->transitions,
	       sizeof(Transition)*i);
	memcpy(&new_transitions[i+1], &trie->transitions[i],
	       sizeof(Transition)*(trie->num_transitions-i));
	free(trie->transitions);
	trie->transitions = new_transitions;
	new_transitions = NULL;
	trie->num_transitions += 1;

	/* Initialize the new transition. */
	transition = &trie->transitions[i];
	transition->suffix = new_suffix;
	transition->next = newtrie;
	transition->next->value = (void *)value;

	if(0) {
	insert_memerror:
	    if(new_transitions) free(new_transitions);
	    if(newtrie) free(newtrie);
	    if(new_suffix) free_suffix(new_suffix);
	    return 1;
	}
    } 
    /* There are three cases where the key and suffix share some
       letters. 
       1.  suffix is proper substring of key.
       2.  key is proper substring of suffix.
       3.  neither is proper substring of other.

       For cases 2 and 3, I need to first split up the transition
       based on the number of characters shared.  Then, I can insert
       the rest of the key into the next trie.
    */
    else {
	/* Count the number of characters shared between key
	   and suffix. */
	int chars_shared = 0;
	while(key[chars_shared] && key[chars_shared] == suffix[chars_shared])
	    chars_shared++;

	/* Case 2 or 3, split this sucker! */
	if(chars_shared < strlen(suffix)) {
	    Trie* newtrie=NULL;
	    char *new_suffix1=NULL, *new_suffix2=NULL;

	    if(!(new_suffix1 = shared_suffix(key, chars_shared)))
		goto split_memerror;
	    if(!(new_suffix2 = duplicate(suffix+chars_shared)))
		goto split_memerror;
	    if(!(newtrie = Trie_new()))
		goto split_memerror;
	    if(!(newtrie->transitions = malloc(sizeof(Transition))))
		goto split_memerror;
	    newtrie->num_transitions = 1;
	    newtrie->transitions[0].next = transition->next;
	    newtrie->transitions[0].suffix = new_suffix2;

	    free_suffix(transition->suffix);
	    transition->suffix = new_suffix1;
	    transition->next = newtrie;

	    if(0) {
	    split_memerror:
		if(newtrie && newtrie->transitions) free(newtrie->transitions);
		if(newtrie) free(newtrie);
		if(new_suffix2) free_suffix(new_suffix2);
		if(new_suffix1) free_suffix(new_suffix1);
		return 1;
	    }
	}
	retval = Trie_set(transition->next, key+chars_shared, value);
    }

    return retval;
}

void Trie_del(Trie* trie) {
    int i;
    if(!trie)
	return;
    for(i=0; i<trie->num_transitions; i++) {
	Transition* transition = &trie->transitions[i];
	if(transition->suffix)
	    free_suffix(transition->suffix);
	Trie_del(transition->next);
    }
    free(trie);
}

void *Trie_get(const Trie* trie, const char *key) {
    int first, last, mid;

    if(!key[0]) {
	return trie->value;
    }

    /* The transitions are stored in alphabetical order.  Do a binary
     * search to find the proper one.
     */
    first = 0;
    last = trie->num_transitions-1;
    while(first <= last) {
	    Transition* transition;
	    char *suffix;
	    int c;
	    mid = (first+last)/2;
	    transition = &trie->transitions[mid];
	    suffix = transition->suffix;
	    /* If suffix is a substring of key, then get the value from
	       the next trie.
	    */
	    c = strncmp(key, suffix, strlen(suffix));
	    if(c < 0)
	        last = mid-1;
	    else if(c > 0)
	        first = mid+1;
	    else
	        return Trie_get(transition->next, key+strlen(suffix));
    }
    return NULL;
}


/* Mutually recursive, so need to make a forward declaration. */
static void
_get_approximate_trie(const Trie* trie, const char *key, const int k,
		      void (*callback)(const char *key, 
				       const void *value,
				       const int mismatches,
				       void *data),
		      void *data, 
		      const int mismatches,
		      char *current_key, const int max_key
		      );

static void 
_get_approximate_transition(const char *key, 
			    const int k,
			    const Transition* transition, 
			    const char *suffix,
			    void (*callback)(const char *key, 
					     const void *value,
					     const int mismatches,
					     void *data),
			    void *data, 
			    const int mismatches,
			    char *current_key, const int max_key
			    )
{
    int i;
    int prev_keylen = strlen(current_key);

    /* Short circuit optimization.  If there's too many characters to
       possibly be a match, then don't even try to match things. */
    if((int)(strlen(suffix) - strlen(key)) > k)
	return;

    /* Match as many characters as possible. */
    i = 0;
    while(suffix[i] && (key[i] == suffix[i])) {
	i++;
    }
    /* Check to make sure the key is not too long.  BUG: If it is,
       fails silently. */
    if((prev_keylen+i) >= max_key)
	return;
    strncat(current_key, suffix, i);

    /* If all the letters in the suffix matched, then move to the
       next trie. */
    if(!suffix[i]) {
	_get_approximate_trie(transition->next, &key[i], k, callback, data,
			      mismatches, current_key, max_key);
    }
    /* Otherwise, try out different kinds of mismatches. */
    else if(k) {
	int new_keylen = prev_keylen+i;

	/* Letter replacement, skip the next letter in both the key and
	   suffix. */
	if((new_keylen+1 < max_key) && key[i] && suffix[i]) {
	    current_key[new_keylen] = suffix[i];
	    current_key[new_keylen+1] = 0;
	    _get_approximate_transition(&key[i+1], k-1,
					transition, &suffix[i+1],
					callback, data,
					mismatches+1, current_key, max_key);
	    current_key[new_keylen] = 0;
	}

	/* Insertion in key, skip the next letter in the key. */
	if(key[i]) {
	    _get_approximate_transition(&key[i+1], k-1, 
					transition, &suffix[i],
					callback, data,
					mismatches+1, current_key, max_key);
	}

	/* Deletion from key, skip the next letter in the suffix. */
	if((new_keylen+1 < max_key) && suffix[i]) {
	    current_key[new_keylen] = suffix[i];
	    current_key[new_keylen+1] = 0;
	    _get_approximate_transition(&key[i], k-1,
					transition, &suffix[i+1],
					callback, data,
					mismatches+1, current_key, max_key);
	    current_key[new_keylen] = 0;
	}
    }
    current_key[prev_keylen] = 0;
}

static void
_get_approximate_trie(const Trie* trie, const char *key, const int k,
		      void (*callback)(const char *key, 
				       const void *value,
				       const int mismatches,
				       void *data),
		      void *data, 
		      const int mismatches,
		      char *current_key, const int max_key
		      )
{
    int i;

    /* If there's no more key to match, then I'm done. */
    if(!key[0]) {
	if(trie->value)
	    (*callback)(current_key, trie->value, mismatches, data);
    }
    /* If there are no more mismatches allowed, then fall back to the
       faster Trie_get. */
    else if(!k) {
	void *value = Trie_get(trie, key);
	if(value) {
	    int l = strlen(current_key);
	    /* Make sure I have enough space for the full key. */
	    if(l + strlen(key) < max_key) {
		strcat(current_key, key);
		(*callback)(current_key, value, mismatches, data);
		current_key[l] = 0;
	    }
	    /* BUG: Ran out of space for the key.  This fails
	       silently, but should signal an error. */
	}
    }
    /* If there are no more transitions, then all the characters left
       in the key are mismatches. */
    else if(!trie->num_transitions) {
	if(trie->value && (strlen(key) <= k)) {
	    (*callback)(current_key, trie->value, 
			mismatches+strlen(key), data);
	}
    }
    /* Otherwise, try to match each of the transitions. */
    else {
	for(i=0; i<trie->num_transitions; i++) {
	    Transition* transition = &trie->transitions[i];
	    const char *suffix = transition->suffix;
	    _get_approximate_transition(key, k, transition, suffix,
					callback, data, 
					mismatches, current_key, max_key);
	}
    }

}


void 
Trie_get_approximate(const Trie* trie, const char *key, const int k,
		     void (*callback)(const char *key, 
				      const void *value,
				      const int mismatches,
				      void *data),
		     void *data
		     )
{
    KEY[0] = 0;
    _get_approximate_trie(trie, key, k, callback, data, 0, KEY,MAX_KEY_LENGTH);
}

int Trie_len(const Trie* trie) 
{
    int i;
    int length = 0;
    
    if(!trie)
	return 0;
    if(trie->value)
	length += 1;
    for(i=0; i<trie->num_transitions; i++) {
	length += Trie_len(trie->transitions[i].next);
    }
    return length;
}

int Trie_has_key(const Trie* trie, const char *key) 
{
    return Trie_get(trie, key) != NULL;
}

int Trie_has_prefix(const Trie* trie, const char *prefix) 
{
    int first, last, mid;

    if(!prefix[0]) {
	return 1;
    }

    /* The transitions are stored in alphabetical order.  Do a binary
     * search to find the proper one.
     */
    first = 0;
    last = trie->num_transitions-1;
    while(first <= last) {
	Transition* transition;
	char *suffix;
	int suffixlen, prefixlen, minlen;
	int c;
	mid = (first+last)/2;
	transition = &trie->transitions[mid];
	suffix = transition->suffix;
	suffixlen = strlen(suffix);
	prefixlen = strlen(prefix);
	minlen = (suffixlen < prefixlen) ? suffixlen : prefixlen;
	c = strncmp(prefix, suffix, minlen);
	if(c < 0)
	    last = mid-1;
	else if(c > 0)
	    first = mid+1;
	else
	    return Trie_has_prefix(transition->next, prefix+minlen);
    }
    return 0;
}

static void 
_iterate_helper(const Trie* trie, 
		void (*callback)(const char *key, 
				 const void *value,
				 void *data),
		void *data,
		char *current_key, const int max_key)
{
    int i;
    if(trie->value)
	(*callback)(current_key, trie->value, data);
    for(i=0; i<trie->num_transitions; i++) {
	Transition* transition = &trie->transitions[i];
	const char *suffix = transition->suffix;
	int keylen = strlen(current_key);

	if(keylen + strlen(suffix) >= max_key) {
	    /* BUG: This will fail silently.  It should raise some
	       sort of error. */
	    continue;
	}
	strcat(current_key, suffix);
	_iterate_helper(transition->next, callback, data, 
			current_key, max_key);
	current_key[keylen] = 0;
    }
}

void 
Trie_iterate(const Trie* trie, 
	     void (*callback)(const char *key, 
			      const void *value,
			      void *data),
	     void *data)
{
    KEY[0] = 0;
    _iterate_helper(trie, callback, data, KEY, MAX_KEY_LENGTH);
}

static void
_with_prefix_helper(const Trie* trie, const char *prefix,
		    void (*callback)(const char *key, 
				     const void *value,
				     void *data),
		    void *data,
		    char *current_key, const int max_key)
{
    int first, last, mid;

    if(!prefix[0]) {
	_iterate_helper(trie, callback, data, current_key, max_key);
	return;
    }

    /* The transitions are stored in alphabetical order.  Do a binary
     * search to find the proper one.
     */
    first = 0;
    last = trie->num_transitions-1;
    while(first <= last) {
	Transition* transition;
	const char *suffix;
	int suffixlen, prefixlen, minlen;
	int c;
	mid = (first+last)/2;
	transition = &trie->transitions[mid];
	suffix = transition->suffix;
	suffixlen = strlen(suffix);
	prefixlen = strlen(prefix);
	minlen = (suffixlen < prefixlen) ? suffixlen : prefixlen;
	c = strncmp(prefix, suffix, minlen);
	if(c < 0)
	    last = mid-1;
	else if(c > 0)
	    first = mid+1;
	else {
	    int keylen = strlen(current_key);
	    if(keylen + minlen >= max_key) {
		/* BUG: This will fail silently.  It should raise some
		   sort of error. */
		break;
	    }
	    strncat(current_key, suffix, minlen);
	    _with_prefix_helper(transition->next, prefix+minlen,
				callback, data, current_key, max_key);
	    current_key[keylen] = 0;
	    break;
	}
    }
}

void 
Trie_with_prefix(const Trie* trie, const char *prefix,
		 void (*callback)(const char *key, 
				  const void *value,
				  void *data),
		 void *data
		 )
{
    KEY[0] = 0;
    _with_prefix_helper(trie, prefix, callback, data, KEY, MAX_KEY_LENGTH);
}



/* Need to declare _serialize_transition here so it can be called from
   _serialize_trie. */
static int _serialize_transition(const Transition* transition, 
			  int (*write)(const void *towrite, const int length,
				       void *data),
			  int (*write_value)(const void *value, void *data),
			  void *data);

/* This library also provides code for flattening tries so that they
 * can be saved and read back in later.  The format of a serialized
 * trie is:
 * TYPE        NBYTES    DESCRIPTION
 * byte        1         Whether or not there is a value
 * variable    variable  If there is a value, let the client store it.
 * byte        1         Number of transitions for this Trie.
 * transition  variable
 *   int       4         Number of characters in the suffix.
 *   suffix    variable  the suffix for this transition
 *   byte      1         Whether or not there is a trie
 *   trie      variable  Recursively points to another trie.
 * 
 * The number of bytes and the endian may vary from platform to
 * platform.
 */

static
int _serialize_trie(const Trie* trie, 
		    int (*write)(const void *towrite, const int length,
				 void *data),
		    int (*write_value)(const void *value, void *data),
		    void *data)
{
    int i;
    unsigned char has_value;

    has_value = (trie->value != NULL);
    if(!(*write)(&has_value, sizeof(has_value), data))
	return 0;
    if(has_value) {
	if(!(*write_value)(trie->value, data))
	    return 0;
    }

    if(!(*write)(&trie->num_transitions, sizeof(trie->num_transitions), data))
	return 0;
    for(i=0; i<trie->num_transitions; i++) {
	if(!_serialize_transition(&trie->transitions[i], 
				  write, write_value, data))
	    return 0;
    }

    return 1;
}

static
int _serialize_transition(const Transition* transition, 
			  int (*write)(const void *towrite, const int length,
				       void *data),
			  int (*write_value)(const void *value, void *data),
			  void *data)
{
    int suffixlen;
    unsigned char has_trie;

    suffixlen = strlen(transition->suffix);
    if(!(*write)(&suffixlen, sizeof(suffixlen), data))
	return 0;
    if(!(*write)(transition->suffix, suffixlen, data))
	return 0;

    has_trie = (transition->next != NULL);
    if(!(*write)(&has_trie, sizeof(has_trie), data))
	return 0;
    if(has_trie) {
	if(!_serialize_trie(transition->next, write, write_value, data))
	    return 0;
    }
    return 1;
}

int Trie_serialize(const Trie* trie, 
		   int (*write)(const void *towrite, const int length, 
				void *data),
		   int (*write_value)(const void *value, void *data),
		   void *data)
{
    int success = _serialize_trie(trie, write, write_value, data);
    (*write)(NULL, 0, data);
    return success;
}

static
int _deserialize_transition(Transition* transition,
			    int (*read)(void *wasread, const int length, 
					void *data),
			    void *(*read_value)(void *data),
			    void *data);

static
int _deserialize_trie(Trie* trie, 
		      int (*read)(void *wasread, const int length, void *data),
		      void *(*read_value)(void *data),
		      void *data)
{
    int i;
    unsigned char has_value;

    if(!(*read)(&has_value, sizeof(has_value), data))
	goto _deserialize_trie_error;
    if(has_value != 0 && has_value != 1)
	goto _deserialize_trie_error;
    if(has_value) {
	if(!(trie->value = (*read_value)(data)))
	    goto _deserialize_trie_error;
    }
    if(!(*read)(&trie->num_transitions, sizeof(trie->num_transitions), data))
	goto _deserialize_trie_error;
    if(!(trie->transitions = 
	 malloc(trie->num_transitions*sizeof(Transition))))
	goto _deserialize_trie_error;
    for(i=0; i<trie->num_transitions; i++) {
	if(!_deserialize_transition(&trie->transitions[i], 
				    read, read_value, data))
	    goto _deserialize_trie_error;
    }
    return 1;
   
 _deserialize_trie_error:
    trie->num_transitions = 0;
    if(trie->transitions) {
	free(trie->transitions);
	trie->transitions = NULL;
    }
    trie->value = NULL;
    return 0;
}

static
int _deserialize_transition(Transition* transition,
			    int (*read)(void *wasread, const int length, 
					void *data),
			    void *(*read_value)(void *data),
			    void *data)
{
    int suffixlen;
    unsigned char has_trie;
    
    if(!(*read)(&suffixlen, sizeof(suffixlen), data))
	goto _deserialize_transition_error;
    if(suffixlen < 0 || suffixlen >= MAX_KEY_LENGTH)
	goto _deserialize_transition_error;
    if(!(*read)(KEY, suffixlen, data))
	goto _deserialize_transition_error;
    KEY[suffixlen] = 0;
    if(!(transition->suffix = duplicate(KEY)))
	goto _deserialize_transition_error;
    if(!(*read)(&has_trie, sizeof(has_trie), data))
	goto _deserialize_transition_error;
    if(has_trie != 0 && has_trie != 1)
	goto _deserialize_transition_error;
    if(has_trie) {
	transition->next = Trie_new();
	if(!_deserialize_trie(transition->next, read, read_value, data))
	    goto _deserialize_transition_error;
    }
    return 1;

 _deserialize_transition_error:
    if(transition->suffix) {
	free_suffix(transition->suffix);
	transition->suffix = NULL;
    }
    if(transition->next) {
	Trie_del(transition->next);
	transition->next = NULL;
    }
    return 0;
}

Trie* Trie_deserialize(int (*read)(void *wasread, const int length, void *data),
		      void *(*read_value)(void *data),
		      void *data)
{
    Trie* trie = Trie_new();
    if(!_deserialize_trie(trie, read, read_value, data)) {
	Trie_del(trie);
	return NULL;
    }
    return trie;
}
