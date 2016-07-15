#include<stdio.h>
#include<stdlib.h>
#include<time.h>

int main()
{
    srand(time(NULL));
    
    int nSamples, nChan, gSize, baseline;
    printf("Enter nSamples, nChan, gSize and baseline to generate random data. It's suggest that baseline >= 0.4 gSize, nSamples >= 0.2*gSize*gSize.\n");
    scanf("%d%d%d%d", &nSamples, &nChan, &gSize, &baseline);
    
    FILE *input_data;
    if ( (input_data=fopen("input.dat","w"))==NULL )
    {
        printf("Cannot open file input.data\n");
        return 0;
    }
    
    fprintf(input_data, "nSamples=%d\n", nSamples);
    fprintf(input_data, "wSize=33\n");
    fprintf(input_data, "nChan=%d\n", nChan);
    fprintf(input_data, "gSize=%d\n", gSize);
    fprintf(input_data, "baseline=%d\n", baseline);
    fprintf(input_data, "cellSize=5.000000\n");
    
    fclose(input_data);
    
    FILE *random_data;
    if ( (random_data=fopen("randnum.dat","wb"))==NULL )
    {
        printf("Cannot open file randnum.data\n");
        return 233;
    }  
    
    int write_num = nSamples * (3 + nChan);
    for (int i = 0; i < write_num; i++)
    {
        int rnd = rand() + 1;
        double data = (double)(rnd) / (double)(RAND_MAX);
        if (fwrite(&data, sizeof(double), 1, random_data) != 1)
        {
            printf("Random data file write error!\n"); 
            return 233;
        }
    }
    
    printf("Random data generation over.");
    return 0;
}