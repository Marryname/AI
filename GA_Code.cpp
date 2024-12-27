/*
�Ŵ��㷨����
*/
//#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*����*/
#define POPSIZE 100       /*��Ⱥ��ģ*/
#define MAXGENS 500     /*����������*/
#define NVARS 2          /*��������*/
#define PCROSSOVER 0.7       /*�ӽ�����*/
#define PMUTATION 0.07   /*�������*/
#define NUM_THREADS 5




int generation;          /*Ŀǰ����*/
FILE *galog;             /*����ļ�*/
FILE *output;                     /* an output(1) file */

struct genotype          /*�������ͣ�GT����Ⱥ��һ����Ա*/
{
	double gene[NVARS];  /*������*/
	double fitness;      /*��Ӧֵ*/
	double rfitness;     /*�����Ӧֵ*/
};

struct boundtype         /*��������*/
{
	double upper;        /*�����Ͻ�*/
	double lower;        /*�����½�*/
};

struct genotype population[POPSIZE + 1];     /*һ����Ⱥ*/
int cur_best;            /*��Ѹ����±�*/
struct genotype newpopulation[POPSIZE + 1];  /*�µ���Ⱥ*/
struct boundtype bounds[NVARS];              /*���ڴ洢ÿ�����������½�*/

/*������������*/
void initialize(void);
double randval(double, double);
void evaluate(void);
void keep_the_best(void);
void elitist(void);
void select(void);
void crossover(void);
void Xover(int, int);
void swap(double*, double*);
void mutate(void);
void report(void);


void main(void) 
{
	int i;
	if ((galog = fopen("galog.txt", "w")) == NULL) 
	{
		exit(1);
	}
	generation = 0;
	fprintf(galog, "\ngeneration |  best  | average | standard\n");
	fprintf(galog, "number     |  value | fitness | deviation\n");

	//omp_set_num_threads(NUM_THREADS);

	initialize();
	evaluate();
	keep_the_best();
	while (generation<MAXGENS)
	{
		generation++;
		select();
		crossover();
		mutate();
		report();
		evaluate();
		elitist();
	}
	fprintf(galog, "\n\n Simulation completed\n");
	fprintf(galog, "\n\n Best member:\n");
	for (i = 0; i < NVARS; i++) 
	{
		fprintf(galog, "\n var(%d)=%3.3f", i, population[POPSIZE].gene[i]);
	}
	fprintf(galog, "\n\n Best fitness =%3.3f", population[POPSIZE].fitness);
	fclose(galog);
	printf("Success\n");
}

/*initialize�����������������ʼ��Ⱥ�塣
��gadata.txt�ж�ȡ���½��ޡ�*/
void initialize(void) 
{
	FILE *infile;
	int i, j;
	double lbound, ubound;

	if ((infile = fopen("gadata.txt", "r")) == NULL) 
	{
		printf( "\nCannot open input file!\n");
		exit(1);
	}
	remove("output.dat");
	for (i = 0; i < NVARS; i++) 
	{
		fscanf(infile, "%lf", &lbound);
		fscanf(infile, "%lf", &ubound);
		/*ÿ���������Լ������½磬�ǿ�����һ��ȫ�ֵ�������߽ṹ��ʾ������ÿ�����嶼�洢�����ӿռ�������*/
		bounds[i].lower = lbound;
		bounds[i].upper = ubound;
//#pragma omp for
		for (j = 0; j < POPSIZE; j++) 
		{
			population[j].fitness = 0;
			population[j].rfitness = 0;
			population[j].gene[i] = randval(bounds[i].lower, bounds[i].upper);
		}
	}
	fclose(infile);
}

/*randval�����ڽ����ڵ�һ��ֵ*/
double randval(double low, double high) 
{
	double val;
	val = ((double)(rand() % 1000) / 1000.0)*(high - low) + low;
	return val;
}

/*evaluate���������û�����
ÿ�θı䶼��Ҫ���±���
������x1^2-x1+x2+x3*/

void evaluate(void) 
{
//	int mem;
//	int i;
//	double x[NVARS + 1];
//#pragma omp for
//	for (mem = 0; mem < POPSIZE; mem++) 
//	{
//		for (i = 0; i < NVARS; i++) 
//			x[i + 1] = population[mem].gene[i];
//		
//		population[mem].fitness = x[1] * x[1] - x[1] * x[2] + x[3];   /*��Ӧֵ�ı�ֻ��Ҫ�޸Ĵ˴�*/
//	}

	if ((output = fopen("output.dat", "a")) == NULL)
	{
		exit(1);
	}
	int mem;
	int i;
	double x[NVARS + 1];

	for (mem = 0; mem < POPSIZE; mem++)
	{
		for (i = 0; i < NVARS; i++)
			x[i + 1] = population[mem].gene[i];

		population[mem].fitness = x[1] * x[1] + x[2] * x[2];


		if (population[mem].fitness <40 && population[mem].fitness>35)
		{
			fprintf(output, "\n%5d,      %6.9f, %6.9f, %6.9f ", generation,
				x[1], x[2], population[mem].fitness);
		}

	}
}

/*keep_the_best��������¼��Ⱥ����Ѹ�����population[POPSIZE]��*/
void keep_the_best() 
{
	int mem;
	int i;
	cur_best = 0;/*�洢���Ÿ�������*/
//#pragma omp for
	for (mem = 0; mem < POPSIZE; mem++) 
	{
		if (population[mem].fitness > population[POPSIZE].fitness) 
		{
			cur_best = mem;
			population[POPSIZE].fitness = population[mem].fitness;
		}
	}
	for (i = 0; i < NVARS; i++)
		population[POPSIZE].gene[i] = population[cur_best].gene[i];
}


/*elitist��������һ������Ѹ���ᱻ�洢�������У�
�����һ������Ѹ������һ������Ѹ���
���߻������ǰ��Ⱥ�е�������*/
void elitist() 
{
	int i;
	double best, worst;
	int best_mem, worst_mem;
	best = population[0].fitness;
	worst = population[0].fitness;
//#pragma omp for
	for (i = 0; i < POPSIZE ; ++i)
	{
		if (population[i].fitness >= best) 
		{
			best = population[i].fitness;
			best_mem = i;

		}
		if (population[i].fitness <= worst) 
		{
			worst = population[i].fitness;
			worst_mem = i;
		}
	}
	if (best >= population[POPSIZE].fitness)
	{
		for (i = 0; i < NVARS; i++)
			population[POPSIZE].gene[i] = population[best_mem].gene[i];
		population[POPSIZE].fitness = population[best_mem].fitness;
	}
	else 
	{
		for (i = 0; i < NVARS; i++)
			population[worst_mem].gene[i] = population[POPSIZE].gene[i];
		population[worst_mem].fitness = population[POPSIZE].fitness;
	}
}

void select(void) 
{
	int mem, i, j, k;
	double sum = 0; 
	double  cfitness;
	double p;
	for (mem = 0; mem < POPSIZE; mem++) 
	    sum += population[mem].fitness;
//#pragma omp for
	for (mem = 0; mem < POPSIZE; mem++)
		population[mem].rfitness = population[mem].fitness / sum;


//#pragma omp for
	for (i = 0; i < POPSIZE; i++) 
	{
		cfitness = 0;
		p = rand() % 1000 / 1000.0;

		for (j = 0; j < POPSIZE; j++)
		{
			cfitness += population[j].rfitness;
			if (p <= cfitness)
			{
				newpopulation[i] = population[j];
				break;
			}
		}
		
	}
//#pragma omp for
	for (i = 0; i < POPSIZE; i++)
		population[i] = newpopulation[i];
}

/*crossover����ѡ�����ӣ�ѡ�������������뽻�棬
ʵ��һ�����㽻��*/
void crossover(void) 
{
	int i, mem, one;
	int first = 0;
	double x;
	for (mem = 0; mem < POPSIZE; ++mem) 
	{
		x = rand() % 1000 / 1000.0;
		if (x < PCROSSOVER) 
		{
			++first;
			if (first % 2 == 0)
				Xover(one, mem);
			else 
				one = mem;
		}
	}
}

/*Crossover:ʵ��������ѡ��ĸ�����Ľ���*/
void Xover(int one, int two) 
{
	int i;
	int point;
	if (NVARS > 1) 
	{
		if (NVARS == 2)
			point = 1;
		else
			point = (rand() % (NVARS - 1)) + 1;/*����2���÷ֳ���������*/
		for (i = 0; i < point; i++)
			swap(&population[one].gene[i], &population[two].gene[i]);
	}
}

/*swap:����������*/
void swap(double*x, double*y) 
{
	double temp;
	temp = *x;
	*x = *y;
	*y = temp;
}

/*Mutation:��������졣
һ��������ѡΪ������һ�������½����ڵ����ֵ�����*/
void mutate(void) 
{
	int i,j;
	double lbound, hbound;
	double x;
//#pragma omp for
	for (i = 0; i < POPSIZE; i++) 
		for (j = 0; j < NVARS; j++) 
		{
			x = rand() % 1000 / 1000.0;
			if (x < PMUTATION) 
			{
				lbound = bounds[j].lower;
				hbound = bounds[j].upper;
				population[i].gene[j] = randval(lbound, hbound);
			}
		}
}
/*report����������ģ��Ĺ��̣��������������ļ��У����ö��š���������*/
void report(void) 
{
	int i;
	double best_val;       /*�����Ⱥ��Ӧֵ*/
	double avg;            /*ƽ����Ⱥ��Ӧֵ*/
	double stddev;         /*��Ⱥ��Ӧֵ�ı�׼ƫ��*/
	double sum_square;     /*��׼�����ƽ����*/
	double square_sum;     /*��׼����ĺ͵�ƽ��*/
	double sum;            /*ȫ����Ⱥ��Ӧֵ*/
	
	sum = 0.0;
	sum_square = 0.0;
	for(i=0;i<POPSIZE;i++)
	{
		sum += population[i].fitness;
		sum_square += population[i].fitness*population[i].fitness;
	}
	avg = sum / (double)POPSIZE;
	square_sum = avg*avg*(double)POPSIZE;
	stddev = sqrt((sum_square - square_sum) / (POPSIZE - 1));
	best_val = population[POPSIZE].fitness;
	fprintf(galog, "\n%7d,   | %6.3f | %6.3f | %6.3f ",
		generation, best_val, avg, stddev);

}

