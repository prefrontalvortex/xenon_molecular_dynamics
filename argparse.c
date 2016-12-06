//
// Created by mm on 12/5/16.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "argparse.h"
#include "aux.h"

void parser_add_arg(arg_t **root, char *argvname, char *type, int position, int ispositional, int isnecessary) {
    arg_t **p_curr, *p_last;
    arg_t *node = (arg_t *) emalloc(sizeof(arg_t));
    node->argvname = (char *) emalloc(strlen(argvname) + 1);
    node->type = (char *) emalloc(strlen(type) + 1);
    node->argvname = argvname;
    node->type = type;
    node->position = position;
    node->ispositional = ispositional;
    node->isnecessary = isnecessary;
    node->prev = NULL;


//    CHECK(9, "start traversal 1");
    p_last = NULL;
    for (p_curr = &(*root); *p_curr != NULL; p_curr = &(*p_curr)->next) {
        if ((*p_curr)->next == NULL) {
            node->prev = (*p_curr);
        }
    } // traverse
//    CHECK(9, "now setting node");
    *p_curr = node;

}

void parser_populate(arg_t **root, int argc, char **argv) {
    int i;
    for (i = 0; i < argc; i++) {
        parser_add_arg(root, argv[i], "arg", i, 0, 0);
    }
}


void ll_traverse(arg_t **root) {
    arg_t **p_curr;

    for (p_curr = &(*root); *p_curr != NULL; p_curr = &(*p_curr)->next) {
        fprintf(stderr, " %s",  (*p_curr)->argvname);
    } // traverse
    fprintf(stderr, "\n");
}


arg_t *ll_index(arg_t **root, int index) {
    arg_t **p_curr;
    int i = 0;
    for (p_curr = &(*root); *p_curr != NULL; p_curr = &(*p_curr)->next) {
//        CHECK(9, "ll index");
        if (i == index) {
            return (*p_curr);
        }
        i++;
    } // traverse
    return NULL;
}

arg_t *ll_search_name(arg_t **root, char *name) {
    arg_t **p_curr;

    for (p_curr = &(*root); *p_curr != NULL; p_curr = &(*p_curr)->next) {
//        CHECK(9, "search");
        if (!strcmp((*p_curr)->argvname, name)) {
            return (*p_curr);
        }
    } // traverse
    return NULL;
}

char *parser_get_tag(arg_t **root, char *tag) {
    arg_t *node = ll_search_name(root, tag);
    if (node == NULL || node->next == NULL) {
        return NULL;
    }
    return node->next->argvname;
}

//void argparser_construct(parsearg_t *parser, int argc, c) {
//
//}

void parse_assign(void *variable, char *argname, enum parsertypes_t type, arg_t *args, char *defaultval) {
    // this only handles tagged variables at the moment
    void *vp;
    int ival;
    double dval;
    char *csval;
    char *tag = parser_get_tag(&args, argname);
    int usedefault = (tag == NULL) ? 1 : 0;
//    if (usedefault) {

//    } else {
    switch (type) {
        case P_INT:
            ival = (usedefault) ? atoi(defaultval) : atoi(tag);
            vp = &ival;
//            *variable = vp;
            break;
        case P_DOUBLE:
            break;
        case P_CHARS:
            break;
        default:
            break;
    }
//    }

}

// todo: WHAT A MESS. Clean this up.

void parse_assign_d(double *variable, char *argname, arg_t *args, char *defaultval) {
    // this only handles tagged variables at the moment
    char *tag = parser_get_tag(&args, argname);
//    int usedefault = (tag == NULL) ? 1 : 0;
    tag = (tag == NULL) ? defaultval : tag;
    *variable = atof(tag);
}

void parse_assign_e(double *variable, char *argname, arg_t *args, char *defaultval) {
    // this only handles tagged variables at the moment
    char *tag = parser_get_tag(&args, argname);
//    int usedefault = (tag == NULL) ? 1 : 0;
    tag = (tag == NULL) ? defaultval : tag;
    sscanf(tag, "%le", variable);
}

void parse_assign_i(int *variable, char *argname, arg_t *args, char *defaultval) {
    // this only handles tagged variables at the moment
    char *tag = parser_get_tag(&args, argname);
//    int usedefault = (tag == NULL) ? 1 : 0;
    tag = (tag == NULL) ? defaultval : tag;
    *variable = atoi(tag);
}

void parse_assign_cs(char **variable, char *argname, arg_t *args, char *defaultval) {
    // this only handles tagged variables at the moment
    char *tag = parser_get_tag(&args, argname);
//    int usedefault = (tag == NULL) ? 1 : 0;
    tag = (tag == NULL) ? defaultval : tag;
    *variable = tag;
}

void parse_assign_b(int *variable, char *argname, arg_t *args, char *defaultval) {
    // this only handles tagged variables at the moment
//    char *tag = parser_get_tag(&args, argname);
    arg_t *node = ll_search_name(&args, argname);
    if (node == NULL) {
        *variable = atoi(defaultval);
    } else {
        *variable = 1;
    }
}


