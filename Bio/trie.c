#include <stdio.h>    /* printf */
#include <stdlib.h>   /* malloc */
#include <string.h>   /* strcmp, strlen */

#include "trie.h"

static char* duplicate(const char* s) {
    /* Don't use strdup, as it's not ANSI C. */
    char* t = malloc((strlen(s)+1)*sizeof(char));
    if (!t) return NULL;
    strcpy(t, s);
    return t;
}


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
	    if(new_suffix) free(new_suffix);
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

	    if(!(new_suffix1 = malloc(chars_shared+1)))
		goto split_memerror;
	    strncpy(new_suffix1, key, chars_shared);
	    new_suffix1[chars_shared] = 0;
	    if(!(new_suffix2 = duplicate(suffix+chars_shared)))
		goto split_memerror;
	    if(!(newtrie = Trie_new()))
		goto split_memerror;
	    if(!(newtrie->transitions = malloc(sizeof(Transition))))
		goto split_memerror;
	    newtrie->num_transitions = 1;
	    newtrie->transitions[0].next = transition->next;
	    newtrie->transitions[0].suffix = new_suffix2;

	    free(transition->suffix);
	    transition->suffix = new_suffix1;
	    transition->next = newtrie;

	    if(0) {
	    split_memerror:
		if(newtrie && newtrie->transitions) free(newtrie->transitions);
		if(newtrie) free(newtrie);
		if(new_suffix2) free(new_suffix2);
		if(new_suffix1) free(new_suffix1);
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
	    free(transition->suffix);
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
	free(transition->suffix);
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
