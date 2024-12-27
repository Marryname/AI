/*
遗传算法代码
*/
//#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*参数*/
#define POPSIZE 100       /*种群规模*/
#define MAXGENS 500     /*最大迭代次数*/
#define NVARS 2          /*变量个数*/
#define PCROSSOVER 0.7       /*杂交概率*/
#define PMUTATION 0.07   /*变异概率*/
#define NUM_THREADS 5




int generation;          /*目前代数*/
FILE *galog;             /*输出文件*/
FILE *output;                     /* an output(1) file */

struct genotype          /*迭代类型（GT）种群的一个成员*/
{
	double gene[NVARS];  /*变量串*/
	double fitness;      /*适应值*/
	double rfitness;     /*相对适应值*/
};

struct boundtype         /*界限类型*/
{
	double upper;        /*变量上界*/
	double lower;        /*变量下界*/
};

struct genotype population[POPSIZE + 1];     /*一个种群*/
int cur_best;            /*最佳个体下标*/
struct genotype newpopulation[POPSIZE + 1];  /*新的种群*/
struct boundtype bounds[NVARS];              /*用于存储每个变量的上下界*/

/*函数过程声明*/
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

/*initialize函数：用随机函数初始化群体。
在gadata.txt中读取上下界限。*/
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
		/*每个变量有自己的上下界，那可以用一个全局的数组或者结构表示，无需每个个体都存储，增加空间冗余性*/
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

/*randval生成在界限内的一个值*/
double randval(double low, double high) 
{
	double val;
	val = ((double)(rand() % 1000) / 1000.0)*(high - low) + low;
	return val;
}

/*evaluate函数，由用户定义
每次改变都需要重新编译
现在是x1^2-x1+x2+x3*/

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
//		population[mem].fitness = x[1] * x[1] - x[1] * x[2] + x[3];   /*适应值改变只需要修改此处*/
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

/*keep_the_best函数：记录种群的最佳个体于population[POPSIZE]。*/
void keep_the_best() 
{
	int mem;
	int i;
	cur_best = 0;/*存储最优个体的序号*/
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


/*elitist函数：上一代中最佳个体会被存储在数组中，
如果这一代的最佳个体比上一代的最佳个体差，
后者会替代当前种群中的最差个体*/
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

/*crossover交叉选择算子：选择两个父代参与交叉，
实现一个单点交叉*/
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

/*Crossover:实现两个被选择的父代间的交换*/
void Xover(int one, int two) 
{
	int i;
	int point;
	if (NVARS > 1) 
	{
		if (NVARS == 2)
			point = 1;
		else
			point = (rand() % (NVARS - 1)) + 1;/*好像2不用分出来的样子*/
		for (i = 0; i < point; i++)
			swap(&population[one].gene[i], &population[two].gene[i]);
	}
}

/*swap:交换两个数*/
void swap(double*x, double*y) 
{
	double temp;
	temp = *x;
	*x = *y;
	*y = temp;
}

/*Mutation:随机化变异。
一个变量被选为变异则被一个在上下界限内的随机值所替代*/
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
/*report函数：报告模拟的过程，数据输出到输出文件中，且用逗号“，”隔开*/
void report(void) 
{
	int i;
	double best_val;       /*最佳种群适应值*/
	double avg;            /*平均种群适应值*/
	double stddev;         /*种群适应值的标准偏差*/
	double sum_square;     /*标准计算的平方和*/
	double square_sum;     /*标准计算的和的平方*/
	double sum;            /*全体种群适应值*/
	
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

