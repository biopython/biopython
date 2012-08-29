%option noyywrap
COMMENT					#.*\n
NAME					_[^ \t\n]+
LOOP					[Ll][Oo][Oo][Pp]_
DATA					[Dd][Aa][Tt][Aa]_[^ \t\n]+
FREE_VALUE				[^ \t\n]+
SINGLE_QUOTE_VALUE		'[^'\n]*'
DOUBLE_QUOTE_VALUE		\"[^"\n]*\"
SEMICOLON_VALUE			^;(.*\n[^;])*.*\n;	

/* 
Number code of the tokens is:

NAME 			(eg. _entity.id)		1
LOOP 			(loop_)					2
DATA			(eg. data_something)	3
SEMICOLON 		(eg. ;value\n;)	 		4	
DOUBLE QUOTE 	(eg. "value")			5
SINGLE QUOTE 	(eg. 'value')			6
FREE 			(eg. value )			7

*/

%%

{COMMENT}					/* ignore */			

{NAME}						{ return 1; }	

{LOOP}						{ return 2; }					

{DATA}						{ return 3;	}

{SEMICOLON_VALUE}			{ return 4; }	

{DOUBLE_QUOTE_VALUE}		{ return 5; }	

{SINGLE_QUOTE_VALUE}		{ return 6; }	

{FREE_VALUE}				{ return 7; }	

[ \t\n]+					/* ignore */

%%					


void mmcif_set_file(FILE *fp)
{
	yyin=fp;
}	

int mmcif_get_token(void)
{
	extern int yylex(void);
	return yylex();
}

char *mmcif_get_string(void)
{
	return yytext;
}	

