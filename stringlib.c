
#include "stringlib.h"

/*   String Related Subroutines */

char *trim ( char *s )
{
	int i = 0;
	int j = strlen ( s ) - 1;
	int k = 0;
	       
	while ( isspace ( s[i] ) && s[i] != '\0' )
		i++;
		 
	while ( isspace ( s[j] ) && j >= 0 )
      		j--;
		   
	while ( i <= j )
       		s[k++] = s[i++];
	     
     	s[k] = '\0';
       
	return s;
}

char *string_firstword ( char *s )
{
	int i = 0;
	int j = strlen ( s ) - 1;
	int k = 0;
	       
	while ( isspace ( s[i] ) && s[i] != '\0' )
		i++;

	j=i;
	while ( s[j] != '\0') {
		if (isspace(s[j])) {
			j--;
			break;
		}
      		j++;
	}
	
	while ( i <= j )
      		s[k++] = s[i++];
	     
     	s[k] = '\0';

	return s;
}

char *string_secondword ( char *s )
{
	int i = 0;
	int j = strlen ( s ) - 1;
	int k = 0;
	     
        // skip beginning white space if any	
	while ( isspace ( s[i] ) && s[i] != '\0' )
		i++;

	// skip first word
	while ( s[i] != '\0') {
		if (isspace(s[i])) 
			break;
      		i++;
	}

        // skip white space 
	while ( isspace ( s[i] ) && s[i] != '\0' )
		i++;

	// find end of second word
	j=i;
	while ( s[j] != '\0') {
		if (isspace(s[j])) {
			j--;
			break;
		}
      		j++;
	}


	// copy second word to new string	   
	while ( i <= j )
       		s[k++] = s[i++];
	     
     	s[k] = '\0';
       
	return s;
}

char *string_thirdword ( char *s )
{
	int i = 0;
	int j = strlen ( s ) - 1;
	int k = 0;
	     
        // skip beginning white space if any	
	while ( isspace ( s[i] ) && s[i] != '\0' )
		i++;

	// skip first word
	while ( s[i] != '\0') {
		if (isspace(s[i])) 
			break;
      		i++;
	}

        // skip white space 
	while ( isspace ( s[i] ) && s[i] != '\0' )
		i++;

	// skip the second word
	while ( s[i] != '\0') {
		if (isspace(s[i])) {
			break;
		}
      		i++;
	}

        // skip white space 
	while ( isspace ( s[i] ) && s[i] != '\0' )
		i++;

	// find the end of the third word
	j=i;
	while ( s[j] != '\0') {
		if (isspace(s[j])) {
			j--;
			break;
		}
      		j++;
	}

	// copy third word to new string	   
	while ( i <= j )
       		s[k++] = s[i++];
	     
     	s[k] = '\0';
       
	return s;
}

char *string_word ( char *s, int wordNumber)
{
	int i = 0;
	int j = strlen ( s ) - 1;
	int k = 0;
	int word;
	     
	// skip white space if any	
	while ( isspace ( s[i] ) && s[i] != '\0' )
		i++;

	for (word=0;word<(wordNumber-1);word++) {

		// skip  word
		while ( s[i] != '\0') {
			if (isspace(s[i])) 
				break;
      			i++;
		}

		// skip white space if any	
		while ( isspace ( s[i] ) && s[i] != '\0' )
			i++;

	}

	// find end target word
	j=i;
	while ( s[j] != '\0') {
		if (isspace(s[j])) {
			j--;
			break;
		}
      		j++;
	}


	// copy second word to new string	   
	while ( i <= j )
       		s[k++] = s[i++];
	     
     	s[k] = '\0';
       
	return s;
}


