/*
 ============================================================================
 Name        : visbility.c
 Author      : Ikaro Silva
 Version     :
 Copyright   : GPL
 Description : Visibility Graph
 ============================================================================
gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"src/visibility.d" -MT"src/visibility.d" -o "src/visibility.o" "../src/visibility.c"
gcc  -o "visibility"  ./src/visibility.o
*/




#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

long input(void);

/* Global variables. */
long N=0;
double *input_data;	/* input data buffer; allocated and filled by input()
 	 	 	 	 	    x and y series will be interleaved starting with x
 */


static char *help_strings[] = {
		"usage: visibility [OPTIONS ...]\n",
		"where OPTIONS may include:",
		" -h               print this usage summary",
		NULL
};

static void help()
{
	int h;

	fprintf(stderr, help_strings[0]);
	for (h = 1; help_strings[h] != NULL; h++)
		fprintf(stderr, "%s\n", help_strings[h]);
	exit(-1);
}



int main(int argc,char* argv[]){

	char ch;
	while ((ch = getopt(argc,argv,"h"))!=EOF)
		switch(ch){
		case 'h':
			help();
			break;
		default:
			fprintf(stderr,"Unknown option for one: '%s'\n",optarg);
			help();
		}

	argc-= optind;
	argv += optind;

	N=input();

	int a=0,b=0,c=0;
	double ya, yc, yb;
	double ta, tc, tb;
	long samples=N/2; //data is interleaved
	int* count=malloc(samples*sizeof(int));
	double* Pk=malloc(samples*sizeof(double));

	//Calculate the visibility count
	int isVisible=1;
	for(a=0;a<N-3;a=a+2){
		*(count+a)=1;
		*(count+a+2)=1;
		ta=*(input_data+a);
		ya=*(input_data+a+1);
		for(b=a+4;b<N-1;b=b+2){
			isVisible=1;
			tb=*(input_data+b);
			yb=*(input_data+b+1);
			for(c=a+2;c<b;c=c+2){
				tc=*(input_data+c);
				yc=*(input_data+c+1);
				if ( yc > (yb + (ya -yb)*(tb-tc)/(tb-ta)) ){
					//Visibility broken, exit loop inner loop only
					isVisible=0;
					break;
				}
			}
			if(isVisible==1){
				//C is visible to A. Add count to both nodes
				*(count+a)= *(count+a) + 1;
				*(count+c)= *(count+c) + 1;
			}
		}
		//Add to the degree distribution P(k)
		*(Pk + *(count+a) ) = *( Pk + *(count+a)) + 1 ;
	}


	//Normalize Pk by the total number of nodes
	for(a=0;a<samples;a++)
		*(Pk + a) = *(Pk + a)/samples ;

	//Output the degree distribution
	for(a=0;a<samples;a++){
		if( *(Pk + a) != 0)
			printf("%u %f\n",a,*(Pk + a));
	}

	return 0;
}

long input()
{
	long maxdat = 0L, npts = 0L;
	double x,y;
	while (scanf("%f %f",&x,&y) == 2) {
		if (npts >= maxdat) {
			double *s;
			maxdat += 50000;	/* allow the input buffer to grow (the increment is arbitrary) */
			if ((s = realloc(input_data, maxdat * sizeof(double))) == NULL) {
				fprintf(stderr,"insufficient memory, exiting program!");
				exit(-1);
			}
			input_data = s;
		}
		input_data[npts] = x;
		npts++;
		input_data[npts] = y;
		npts++;
	}
	if (npts < 1){
		printf(stderr,"Error, no data read!");
		exit(-1);
	}
	return (npts);
}
