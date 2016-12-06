//
// Created by mm on 12/5/16.
//

#ifndef MOL_ARGPARSE_H
#define MOL_ARGPARSE_H


enum parsertypes_t {
    P_VOID = 0,
    P_INT,
    P_DOUBLE,
    P_CHARS
};

typedef struct _arg_t {
    char *argvname;
    char *type;
    int position;
    int ispositional;
    int isnecessary;
    char *argv;
    struct _arg_t *next;
    struct _arg_t *prev;
} arg_t;

typedef struct _parsearg_t {
    arg_t **input_args;
    arg_t **set_args;
} parsearg_t;

void parser_add_arg(arg_t **root, char *argvname, char *type, int position, int ispositional, int isnecessary);
void parser_populate(arg_t **root, int argc, char **argv);

void ll_traverse(arg_t **root);
arg_t *ll_index(arg_t **root, int index);
arg_t *ll_search_name(arg_t **root, char *name);
char *parser_get_tag(arg_t **root, char *tag);

void parse_assign_d(double *variable, char *argname, arg_t *args, char *defaultval);
void parse_assign_e(double *variable, char *argname, arg_t *args, char *defaultval);
void parse_assign_i(int *variable, char *argname, arg_t *args, char *defaultval);
void parse_assign_cs(char **variable, char *argname, arg_t *args, char *defaultval);
void parse_assign_b(int *variable, char *argname, arg_t *args, char *defaultval);



#endif //MOL_ARGPARSE_H
