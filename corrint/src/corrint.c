/*
 ============================================================================
 Name        : corrint.c
 Author      : Ikaro Silva
 Version     :
 Copyright   : GPL
 Description : Analysis of Time Series based on Correlation Integral
 ============================================================================
 To build:
 gcc  -o corrint  corrint.c   -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

/* Function prototypes. */
long input(void);
void get_err(int windowN, int stepSize,int timeLag, double* err, int nFlag, double th);
void countNeighbors(double *th,double *count, int countN, int Nerr, double* err, int nFlag);
void predictHalf(int windowN, int stepSize,int timeLag, int neighbors, int linearStateStim, int dimSize);
void smooth(int windowN, int stepSize,int timeLag, int neighbors);
void xcorr();
void linearFit(int* blockIndex,int dimSize,int stepSize, int neighboors, double* m,double* b);
/* End of Function prototypes. */


/* Global variables. */
double *input_data;	/* input data buffer; allocated and filled by input() */
long N=0;
long maxLag=0;
double *Rxx; /*autocorrelation of the input data, with length maxLag */
/* End of Global variables. */


static char *help_strings[] = {
		"usage: corrint [OPTIONS ...]\n",
		"where OPTIONS may include:",
		" -h               print this usage summary",
		" -R               calculates autocorrelation only and exit",
		" -p               generate recurrence data only and exit",
		" -d int           embedded dimension size ",
		" -t int           time lag between states (if -1, estimate timeLag from first zero crossing of autocorrelation)",
		" -s int           time lag within state samples",
		" -r double        distance threshold (if 0, default is set to variance of series divided by 10)",
		" -D               Debug Flag, if true prints program detail",
		" -N               Normalize Flag, if true normalize count",
		" -v               Estimates correlation dimension and scaling region given an embedded dimension (-d) parameter",
		" -a               Use this option along with the '-v' option to estimate correlation dimension based on only 2 threshold"
		"                  values determined from the data series\n"
		" -n int           Number of closest neighbors used for prediction ",
		" -P               User first half of the time series as a model to predict the second half (point by point)",
		" -S               Filter the time series by attempting to predict current point  based on the other points",
		" -A               Use mean and slope (linear state approximation) when doing prediction and smoothing (optinos "
		"                  -P and -S.",
		NULL
};

static void help()
{
	int i;
	for (i = 0; help_strings[i] != NULL; i++)
		fprintf(stderr, "%s\n", help_strings[i]);
	exit(-1);
}

int main(int argc,char* argv[]) {

	//Define the parameters of the correlation integral
	int timeLag=2; //Offset with respect only to the first state
	int dim=2;  //Embedded dimension. Limits maximum estimated dimension : v < 2*dim+1
	char ch;
	int stepSize=1;
	int normalizeFlag=1, corrFlag=0, recurFlag=0, estimateDim=0, predictSecondHalf=0, smoothFlag=0;
	int autoEstimateDim=0, linearStateStim=0;
	int windowN;
	register int i;
	//th_arr should be sorted for speed efficiency
	const int countN=19;
	double TH=-1;
	double TH_ARR[19]={0.02, 0.1, 0.2, 0.3, 0.4, 0.5,0.75, 1.25, 1.5, 1.75, 2.0 , 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0};
	double count[countN];
	int neighbors=10;

	for(i=0;i<countN;i++)
		count[i]=0;

	while ((ch = getopt(argc,argv,"hvd:t:s:Rpr:Nn:PSaA"))!=EOF )
		switch(ch){
		case 'v':
			estimateDim=1;
			break;
		case 'n':
			neighbors=atoi(optarg);
			break;
		case 'd':
			dim=atoi(optarg);;
			break;
		case 't':
			timeLag=atoi(optarg);
			break;
		case 's':
			stepSize=atoi(optarg);
			break;
		case 'r':
			TH=atof(optarg);
			break;
		case 'N':
			normalizeFlag=1;
			break;
		case 'R':
			/*Calculate autocorrelation only and exits */
			corrFlag=1;
			break;
		case 'p':
			recurFlag=1;
			break;
		case 'P':
			predictSecondHalf=1;
			break;
		case 'S':
			smoothFlag=1;
			break;
		case 'a':
			autoEstimateDim=1;
			break;
		case 'A':
			linearStateStim=1;
			break;
		case 'h':
			help();
			break;
		default:
			fprintf(stderr,"Unknown option: '%s'\n",optarg);
			help();
			break;
		}

	//Load data into input_data and get the number of samples read
	N=input();

	//Calculate window size here in case dim is entered by user
	windowN=dim*stepSize;

	/* Define single threshold value if passed by user, otherwise use standard array */
	double* th_arr=TH_ARR;
	if(TH != -1 ){
		if(TH != 0 ){
		*th_arr=TH;
		}else{
			//Use variance of the signal to estimate threshold
			double var=0, mx=0;
			for(i=0;i<N;i++){
				mx+= input_data[i];
				var+= input_data[i] * input_data[i];
			}
			mx=mx/((double) N);
			var= var/((double) N) - mx*mx;
			*th_arr= var/10.0;
		}
	}

	if(corrFlag){
		/*Calculate autocorrelation and exit */
		xcorr();
		for(i=0;i<maxLag;i++)
			fprintf(stdout,"%f\n",Rxx[i]);
		exit(0);
	}

	//If we are estimating the signal dimension, we need to set the timeLag to the first zero crossing of the autocorrelation function
	//Otherwise it will be either the default or what the user enters
	if(estimateDim || (timeLag<0)){
		xcorr();
		for(i=0;i<maxLag;i++){
			if(Rxx[i] <=0 ){
				timeLag=i;
				break;
			}
		}
	}
	if(estimateDim){
		//Give a warning early on if we might not have enough point to reliably
		//estimate the dimension given the embedding parameters
		//This is a conservative minimum (sufficient but not necessary).
		double minPoints= pow(10.0,dim);
		if(N< minPoints)
			fprintf(stderr,"Possibly not have enough points to estimate dimension. Total points: %u, sufficient minimum required: %u\n",N,(long) minPoints);

		//Overwrite other parameters accordingly
		normalizeFlag=1;
		recurFlag=0;

	}

	//Find how many distance points we need to calculate and allocate memory accordingly
	int k, errN=0,test=0;
	//TODO: Find  a way to put this in close form.
	// According to Kaplan and Glass this should be:
	//  (N-(1-dim)*stepSize)(N-(1-dim)*stepSize - 1) - timeLag
	//But we should check this properly

	for(i=N-1;i>=windowN-1;i--,test++)
		for(k=i-(windowN-1)-timeLag;k>=windowN-1;k--,errN++);

	double *err;
	if ((err = malloc(errN * sizeof(double))) == NULL) {
		fprintf(stderr,"corrint: insufficient memory for error matrix, exiting program\n");
		exit(-1);
	}

	//Calculate the error matrix
	get_err(windowN,stepSize,timeLag,err,recurFlag,*th_arr);

	//If in prediction mode, predict based or error matrix and exit from here
	if(predictSecondHalf){
		predictHalf(windowN,stepSize,timeLag,neighbors,linearStateStim,dim);
		exit(0);
	}

	//If in smooth mode, smooth based or error matrix and exit from here
	if(smoothFlag){
		smooth(windowN,stepSize,timeLag,neighbors);
		exit(0);
	}

	//Get neighborhood count
	countNeighbors(th_arr,count,countN,errN,err,normalizeFlag);

	//Display results in column format
	if(!estimateDim){
		for(i=0;i<countN;i++){
			fprintf(stdout,"%f \t %f\n",*(th_arr+i),count[i]);
		}
	}

	//For the dimension estimation case, we also print out the estimate slope value
	double slope=0, slope2=0;
	if(estimateDim){
		if(autoEstimateDim != 1){
			slope=(log(count[1]) - log(count[0]) )/( log(*(th_arr+1)) -log(*th_arr) );
			fprintf(stdout,"%f \t %f\n",*th_arr,slope);
			for(i=1;i<countN-1;i++){
				slope=(log(count[i]) - log(count[i-1]) )/( log(*(th_arr+i)) -log(*(th_arr+i-1)) );
				slope2=(log(count[i+1]) - log(count[i]) )/( log(*(th_arr+i+1)) -log(*(th_arr+i)) );
				fprintf(stdout,"%f \t %f\n",*(th_arr+i),(slope+slope2)/2.0);
			}
			slope=(log(count[countN-1]) - log(count[countN-2]) )/( log(*(th_arr+countN-1)) -log(*(th_arr+countN-2)) );
			fprintf(stdout,"%f \t %f\n",*(th_arr+countN-1),slope);
		}else{
			//For the autoEstimateDim case, the first threshold is r1=signma/4 and the second one is
			//the one that yield a count closest to 5x the first
			//Get variance of the time series
			double* r1=malloc(sizeof(double));
			double* countR1=malloc(sizeof(double));
			double var=0, dc=0;
			for(i=0;i<N;i++){
				var+= input_data[i]*input_data[i];
				dc += input_data[i];
			}
			var= var/((double) N) - dc*dc/( (double) N*N);
			*r1=sqrt(var)/4.0;
			countNeighbors(r1,countR1,1,errN,err,normalizeFlag);

			//Get value that is closest to 5x the r1
			double opt=(*countR1)/5.0, val, best=INFINITY;
			int bestInd=-1;
			for(i=0;i<countN;i++){
				val=abs(opt-count[i]);
				if(val<best){
					best=val;
					bestInd=i;
				}
			}
			if(bestInd <0){
				fprintf(stderr,"Could not find optimal scaling region for series ( std=%f, R1=%f, C(R2)=%f, log(C(R2))=%f ).\n",var,*r1,opt,log(opt));
				exit(1);
			}
			slope=(log(*countR1) - log(count[bestInd]) )/( log(*r1) -log(*(th_arr+bestInd)) );
			fprintf(stdout,"%f \t %f\n",*r1,slope);
		}

		//Print estimated lag
		fprintf(stdout,"lag=%u\n",timeLag);
	}

	//Free memory allocated by input
	free(input_data);
	return EXIT_SUCCESS;
}

//Estimate the distance between the states
void get_err(int windowN, int stepSize,int timeLag, double* err, int recurFlag, double th){
	int i, k, z,  index=0;
	double tmpErr;
	//Loop through the data array starting from the end, going to the beginning.
	//For each loop a state vector of size dim and offset timeLag is generated

	for(i=N-1;i>=windowN-1;i--){
		for(k=i-(windowN-1)-timeLag;k>=windowN-1;k--){
			tmpErr=0;
			for(z=0;z<windowN;z+=stepSize){
				tmpErr+=fabs(input_data[i-z]-input_data[k-z]);
			}
			err[index]=tmpErr;

			//for recurrence plots, print any state that has at least one neighbor and exit
			if(recurFlag == 1){
				if(err[index]<th){
					fprintf(stdout,"%u\t%u\n",i,k);
					fprintf(stdout,"%u\t%u\n",k,i);
				}
			}
			index++;
		}
	}
	if(recurFlag){
		exit(0);
	}
}

void smooth(int windowN, int stepSize,int timeLag, int neighbors){
	int i, k, z, n;
	double dist;

	//The iterative approach should be the same as get_err, logging the closest neighboring values and their predictions
	double* blockDistance=calloc(neighbors,sizeof(double));
	int* blockIndex=calloc(neighbors,sizeof(int));
	double maxBlockDistance=-1;;
	int count=0, maxInd;
	double prediction;
	double err=0, cov;
	double dc=0;

	for(i=0;i<N-1;i++){
		dc+=input_data[i];
	}
	dc= dc/( (double) N);


	for(i=windowN;i<N-1;i++){
		//Reset predictions and neighborhood parameters
		for(n=0;n<neighbors;n++){
			*(blockDistance+n)=-1;
			*(blockIndex+n)=0;
		}
		maxBlockDistance=-1;
		maxInd=0;
		count=0;

		for(k=windowN;k<N-1;k++){
			//Get distance from current point if it lies outside its region
			if( k >= i-windowN && k < i)
				continue;
			if( (k-windowN) >= (i-windowN) && (k-windowN) < i)
				continue;
			if( k==i )
				continue;
			dist=0;
			for(z=0;z<windowN;z+=stepSize){
				dist+=fabs(input_data[i-z]-input_data[k-z]);
			}

			//Add point if distance is small or neighborhood is not full
			if(count<neighbors){
				//Filling up the hood
				*(blockDistance+count)=dist;
				*(blockIndex+count)=k;
				if(maxBlockDistance < dist ){
					maxBlockDistance=dist;
					maxInd=count;
				}
				count++;
			}else{
				//If point is cool enough, kick the lamest brother out of the hood
				if(dist<maxBlockDistance){
					*(blockDistance+maxInd)=dist;
					*(blockIndex+maxInd)=k;
					maxBlockDistance=dist;
					//Recalculate the newest lamest bro
					for(n=0;n<neighbors;n++){
						if(maxBlockDistance < *(blockDistance+n) ){
							maxBlockDistance= *(blockDistance+n);
							maxInd=n;
						}
					}
				}
			}
		} //End of neighbor search

		//End of pass, average predictions of all the hood
		prediction=0;
		for(n=0;n<count;n++){
			prediction += ( input_data[ *(blockIndex+n) + 1 ])/count;
		}
		//Output Prediction and true value
		fprintf(stdout,"%f\t%f\n",prediction,input_data[i+1]);

		//Calculate cumulative error and covariance of time series
		err+=(prediction-input_data[i+1])*(prediction-input_data[i+1]);
		cov+=(dc-input_data[i+1])*(dc-input_data[i+1]);
	}

	fprintf(stdout,"err/cov = %f\n",(err/cov));
	exit(0);
}

void predictHalf(int windowN, int stepSize,int timeLag, int neighbors,int linearStateStim, int dimSize){
	int i, k, z, n;
	double dist;

	//The iterative approach should be the same as get_err, logging the closest neighboring values and their predictions
	double* blockDistance=calloc(neighbors,sizeof(double));
	int* blockIndex=calloc(neighbors,sizeof(int));

	double maxBlockDistance=-1;;
	int count=0, maxInd;
	double prediction;
	double err=0, cov;
	double dc=0;

	for(i=N/2+timeLag;i<N-1;i++){
		dc+=input_data[i];
	}
	dc= dc/( (N-1) - (N/2) );


	for(i=N/2+timeLag;i<N-1;i++){
		//Reset predictions and neighborhood parameters
		for(n=0;n<neighbors;n++){
			*(blockDistance+n)=-1;
			*(blockIndex+n)=0;
		}
		maxBlockDistance=-1;
		maxInd=0;
		count=0;

		for(k=(N/2)-timeLag;k>=windowN-1;k--){
			//Get distance from current point
			dist=0;
			for(z=0;z<windowN;z+=stepSize){
				dist+=fabs(input_data[i-z]-input_data[k-z]);
			}

			//Add point if distance is small or neighborhood is not full
			if(count<neighbors){
				//Filling up the hood
				*(blockDistance+count)=dist;
				*(blockIndex+count)=k;
				if(maxBlockDistance < dist ){
					maxBlockDistance=dist;
					maxInd=count;
				}
				count++;
			}else{
				//If point is cool enough, kick the lamest brother out of the hood
				if(dist<maxBlockDistance){
					*(blockDistance+maxInd)=dist;
					*(blockIndex+maxInd)=k;
					maxBlockDistance=dist;
					//Recalculate the newest lamest bro
					for(n=0;n<neighbors;n++){
						if(maxBlockDistance < *(blockDistance+n) ){
							maxBlockDistance= *(blockDistance+n);
							maxInd=n;
						}
					}
				}
			}
		}

		//End of pass, average predictions of all the hood
		prediction=0;
		if(linearStateStim==1){
			double *m,*b;
			b=calloc(dimSize,sizeof(double));
			m=calloc(1,sizeof(double));
			linearFit(blockIndex,dimSize,stepSize,count,m,b);
			prediction= (*m);
			for(n=0;n<dimSize;n++)
				prediction += (*b) * input_data[ i - n*stepSize];
			free(b);
		}else{
			for(n=0;n<count;n++){
				prediction += input_data[ *(blockIndex+n) +1 ]/count;
			}
		}
		//Output Prediction and true value
		fprintf(stdout,"%f\t%f\n",prediction,input_data[i+1]);

		//Calculate cumulative error and covariance of time series
		err+=(prediction-input_data[i+1])*(prediction-input_data[i+1]);
		cov+=(dc-input_data[i+1])*(dc-input_data[i+1]);
	}

	fprintf(stdout,"err/cov = %f\n",(err/cov));
	exit(0);
}

//Get the number of states within a minimum threshold
void countNeighbors(double *th,double *count_arr, int countN, int Nerr, double* err,int nFlag){
	register int i,k;
	//Loop through the distance matrix and then go over the th_arr
	//in decreasing order for each element that is below the threshold,
	//incrementing to the count

	//This assumes that th is sorted in *increasing* size of threshold values
	for(i=0;i<Nerr;i++){
		for(k=0;k<countN;k++){
			if(err[i]<th[countN-1-k]){
				count_arr[countN-1-k]++;
			}else{
				//Lowest possible threshold reached, break from th_arr count
				break;
			}
		}
	}

	//Normalize the count
	if(nFlag==1){
		for(k=0;k<countN;k++)
			count_arr[k]=count_arr[k]/(double) Nerr;
	}

}


/* Estimate autocorrelation */
void xcorr(){
	volatile long lag, n;
	double R, R0;
	maxLag=N/2;
	Rxx=malloc(maxLag * sizeof(double));

	/*Subtract mean from the time series */
	double dc=0;
	for(n=0;n<N;n++)
		dc+=input_data[n];
	dc=(double) dc/N;

	for(n=0;n<N;n++)
		input_data[n]=input_data[n]-dc;

	for(lag=0;lag<maxLag;lag++){
		for(n=0,R=0,R0=0;n<N-maxLag;n++){
			R+= input_data[n+lag]*input_data[n];
			R0+=input_data[n]*input_data[n];
		}
		Rxx[lag]=R/R0;
	}
}


/* Read input data, allocating and filling input_data[].
   The return value is the number of points read.

   This function allows the input buffer to grow as large as necessary, up to
   the available memory (assuming that a long int is large enough to address
   any memory location).
 */
long input()
{
	long maxdat = 0L, npts = 0L;
	double y;
	while (scanf("%lf", &y) == 1) {
		if (npts >= maxdat) {
			double *s;
			maxdat += 50000;	/* allow the input buffer to grow (the increment is arbitrary) */
			if ((s = realloc(input_data, maxdat * sizeof(double))) == NULL) {
				fprintf(stderr,"corrint: insufficient memory, exiting program!");
				exit(-1);
			}
			input_data = s;
		}
		input_data[npts] = y;
		npts++;
	}
	if (npts < 1){
		printf(stderr,"Error, no data read!");
		exit(-1);
	}
	return (npts);
}

void linearFit(int* indeces,int dimSize,int stepSize, int N, double* M, double* b){
	//Do a ordinary least square regression
	int i,k;
	double XX, x;
	*M=0;
	//Calculate the mean and subtract it when calculating the covariances
	for(i=0;i<N;i++)
		(*M)= (*M) + input_data[ indeces[i] + 1]/((double) N);

	//Calculat the coefficient for each variable
	for(i=0;i<dimSize;i++){
		XX=0;
		for(k=0;k<N;k++){
			x=input_data[ indeces[i] - i*stepSize];
			XX += x*x;
			b[i]= x*( input_data[ indeces[i] +1] - (*M) );
		}
		if(XX  == 0){
			b[i]=0;
		}else{
			b[i] = b[i]/XX;
		}
	}
}
