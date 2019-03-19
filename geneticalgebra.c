/********************************************************************
遺伝的アルゴリズム 
********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>

#define POP_SIZE 5   /* 個体数 (必ず奇数に設定) */
#define G_LENGTH 10  /* 個体の遺伝子型のビット数 */
#define MAX_GEN  20  /* 世代数 */
#define M_RATE   0.1 /* 突然変異率 (0?1) */ 
#define MAX_Fitt 10 /* 最高適合率*/ 

/********************************************************************
乱数の発生 (Seedの決定)
********************************************************************/
void init_rnd()
{
	srand((unsigned int)time(NULL));
}

/********************************************************************
乱数の発生 (0?1の乱数)
********************************************************************/
double Random()
{
	return((double)rand() / RAND_MAX);
}

/********************************************************************
遺伝子の初期化
引数 gene[p][i] : 遺伝子pのi番目の成分
********************************************************************/
void init_gene(int gene[POP_SIZE][G_LENGTH])
{
	int p, i;

	/* 乱数の初期化 */
	init_rnd();

	/* 遺伝子を初期化  0?1の乱数を発生し、0.5以上なら1
	0.5未満なら0 */
	printf("<< 初期個体群 >>\n");

	for (i = 0; i < POP_SIZE; i++) {//乱数の追加と仕訳
		for (p = 0; p < G_LENGTH; p++) {
			gene[i][p] = rand()%2;
			printf("%d%d\n", gene[i][p]);
		}
	}

	for (i = 0; i < POP_SIZE; i++) {//数値の整理
		for (p = 0; p < G_LENGTH; p++) {
			if (gene[i][p] < 0.5) {
				gene[i][p] = 0;
			}
			else {
				gene[i][p] = 1;
			}
		}
	}
}

/********************************************************************
適応度の計算
引数 gene[p][i] : 遺伝子pのi番目の成分
fitness[p] : 遺伝子pの適応度
********************************************************************/
void calc_fitness(int gene[POP_SIZE][G_LENGTH], double fitness[POP_SIZE])
{
	int p, i;

	/* 適応度の計算 前半の5bitは0の数 後半の5bitは1の数 */

	for (i = 0; i < POP_SIZE; i++) {
		fitness[i] = 0;//初期成績の追加
		for (p = 0; p < G_LENGTH; p++) {
			if (p < 5 && gene[i][p] == 0) {//適応度の計算
				fitness[i]++;
			}
			else if (p <= 10 && gene[i][p] == 1) {
				fitness[i]++;
			}
		}
	}
}

/**********************************************************************
遺伝子の表示 & 最大適応度・平均適応度の計算 & ファイルへの書き出し
引数 t          : 世代数
gene[p][i] : 遺伝子pのi番目の成分
fitness[p] : 遺伝子pの適応度
*fp        : ファイルポインタ
**********************************************************************/
void show_gene(int t, int gene[POP_SIZE][G_LENGTH],
	double fitness[POP_SIZE],
	FILE *fp)
{
	int p, i;
	double avg_fit=0;  /* 平均適応度  */
	double max_fit=0;  /* 最大適応度  */

					 /* 個体の値、適応度の表示 */
					 /* 平均・最大適応度の計算 */

	printf("適応度表示");
	for (i = 0; i < POP_SIZE; i++) {
		for (p = 0; p < G_LENGTH; p++) {
			printf("%d\n", fitness[i]);
			if (max_fit<fitness[i]) {
				max_fit = gene[i][p];//最大適応の個体群を保存
			}
			avg_fit += gene[i][p];//合算
		}
	}
	avg_fit += avg_fit / POP_SIZE;//平均値の計算
	
	
					 /* 平均・最大適応度の表示 */
	printf("平均適応度 : %lf\n", avg_fit);
	printf("最大適応度 : %lf\n", max_fit);

	/* 平均・最大適応度をファイルに書き込む */
	fprintf(fp, "%d %lf %lf\n", t, avg_fit, max_fit);
}

/**********************************************************************
個体番号 p1 と p2 の適応度と遺伝子を交換
引数 p1, p2     : 遺伝子の番号
gene[p][i] : 遺伝子pのi番目の成分
fitness[p] : 遺伝子pの適応度
**********************************************************************/
void swap_gene(int p1, int p2, int gene[POP_SIZE][G_LENGTH],
	double fitness[POP_SIZE])
{
	int tmp;
	double f;
	int i,p;

	/* 遺伝子型の交換 (遺伝子p1と遺伝子p2の値を入れ替える) */

	while (true) {
		//その世代の個体群の全ての個体の成績（選んだ品物の価値の合計）を算出する
		for (i = 0; i < POP_SIZE; i++) {
			fitness[i] = 0;//初期個体の生成
			for (p = 0; p < G_LENGTH; p++) {
				fitness[i] += gene[i][p];
			}
			if (fitness[i] > MAX_GEN) {
				fitness[i] = 0.0;
			}
		}
	}

	/* 適応度の交換 (遺伝子p1と遺伝子p2の適応度の値を入れ替える) */
}

/**********************************************************************
個体番号 p1 の適応度と遺伝子型を p2 にコピー
引数 p1, p2     : 遺伝子の番号
gene[p][i] : 遺伝子pのi番目の成分
fitness[p] : 遺伝子pの適応度
**********************************************************************/
void copy_gene(int p1, int p2, int gene[POP_SIZE][G_LENGTH],
	double fitness[POP_SIZE])
{
	int i;

	/* 遺伝子のコピー (遺伝子p1を遺伝子p2にコピーする) */
	/* 適応度のコピー (遺伝子p1の適応度を遺伝子p2の適応度にコピーする)*/
	while (true) {
		//その世代の個体群の全ての個体の成績（選んだ品物の価値の合計）を算出する
		for (i = 0; i < POP_SIZE; i++) {
			fitness[i] = 0;//初期個体の生成
			for (p = 0; p < G_LENGTH; p++) {
				fitness[i] +=gene[i][p];
			}
			if (fitness[i] > MAX_GEN) {
				fitness[i] = 0.0;
			}
		}
	}
}

/**********************************************************************
エリート保存
(最小適応度の個体に最大適応度の個体のデータをコピー)
引数 gene[p][i] : 遺伝子pのi番目の成分
fitness[p] : 遺伝子pの適応度
**********************************************************************/
void elite(int gene[POP_SIZE][G_LENGTH], double fitness[POP_SIZE])
{
	int p, i;
	double max_fitness = fitness[0];
	double min_fitness = fitness[0];
	int max_p = 0;
	int min_p = 0;
	int M_F

	/* 最大適応度の個体(max_p)と最小適応度の個体(min_p)を見つける */

	//最大適応度
	for (i = 0; i < POP_SIZE; i++) {
		if (max_fitness < fitness[i]) {
			max_fitness < fitness[i];//成績をソートして最大値を保存
			M_F = i;//エリート選択された番号を記録
		}
	}
	//最小適応度
	for (i = 0; i <POP_SIZE; i++) {
			if (max_fitness < fitness[i]) {//最優秀遺伝子の保存
				min_p = gene[M_F][i];
			}
	}


	/* 最小適応度の個体に最大適応度の個体をコピー */
	copy_gene(max_p, min_p, gene, fitness);
	/* 最大適応度の個体を0番目に移動 */
	swap_gene(0, max_p, gene, fitness);
}

/**********************************************************************
ルーレット選択
引数 gene[p][i] : 遺伝子pのi番目の成分
fitness[p] : 遺伝子pの適応度
**********************************************************************/
void reproduction(int gene[POP_SIZE][G_LENGTH], double fitness[POP_SIZE])
{
	double sum_of_fitness; /* 個体の適応度の総和 */
	double border;         /* ルーレット上の個体間の境界 */
	double r;              /* ルーレット上の選択位置 */
	int p, i;               /* 選ばれた個体の番号 */
	int num;               /* 0 <= num <= POP_SIZE-1 */
	int new_gene[POP_SIZE][G_LENGTH];

	/* ルーレットの1周分 sum_of_fitness を求める */
	sum_of_fitness = 0.0;
	for (p = 0; p<POP_SIZE; p++) {
		sum_of_fitness += fitness[p];
	}

	/* ルーレットを POP_SIZE 回だけ回して次世代の個体を選ぶ */
	for (p = 1; p<POP_SIZE; p++) {
		/* ルーレットを回して場所を選ぶ
		r : 選ばれた位置 (0 <= r <= sum_of_fitness) */
		r = sum_of_fitness * Random();
		/* 選ばれた場所に該当する個体が何番か調べる
		num : 選ばれた個体の番号 (0 <= num <= POP_SIZE-1) */
		num = 0;
		border = fitness[0]; /* 個体間の境界 */
		while (border<r) {
			num++;
			border += fitness[num];
		}

		/* 遺伝子の代入 */
		for (i = 0; i<G_LENGTH; i++) {
			new_gene[p][i] = gene[num][i];
		}
	}

	/* 遺伝子のコピー */
	for (p = 1; p<POP_SIZE; p++) {
		for (i = 0; i<G_LENGTH; i++) {
			gene[p][i] = new_gene[p][i];
		}
	}
}

/**********************************************************************
一点交叉
引数 gene[p][i] : 遺伝子pのi番目の成分
**********************************************************************/
void crossover(int gene[POP_SIZE][G_LENGTH])
{
	int gene1[G_LENGTH]; /* 親1の遺伝子型 */
	int gene2[G_LENGTH]; /* 親2の遺伝子型 */
	int i, j;
	int c_pos;   /* 交叉位置 (1 <= c_pos <= G_LENGTH-1) */

				 /* 交叉位置を1?G_LENGTH-1の範囲でランダムに決め、
				 それより後ろを入れ替える。
				 gene[1]とgene[2],  gene[3]とgene[4] ... のように親にする */

	for (i = 0; i < POP_SIZE;i+=2) {
		c_pos = rand()%G_LENGTH;//切除点をランダムに指定
		for (j = 0; j < G_LENGTH; j++) {
			gene2[i] = gene[i][j];
			gene[i][j] = gene[i + 1][j];//一度次の値に遺伝子を逃がしてから再代入することで交差を実現
			gene[i + 1][j] = gene2[i];
		}
	}

}

/**********************************************************************
突然変異
引数 gene[p][i] : 遺伝子pのi番目の成分
**********************************************************************/
void mutation(int gene[POP_SIZE][G_LENGTH])
{
	int p, i;
	double random;

	/* 0?1の乱数を発生させ、その値が M_RATE 以下ならば
	遺伝子の値をランダムに変える (0ならば1、1ならば0) */

	if (i = 0; i < POP_SIZE; i++) {
		for (p = 0; p < G_LENGTH; p++) {
			random = (double)rand() / 101 / 100;//突然変異率
			if (random < M_RATE) {
				gene[i][p] = rand() % 2;//ランダムな1と0を代入
			}
		}
	}

}

/**********************************************************************
メインプログラム
**********************************************************************/
int main(int argc, char *argv[])
{
	int gene[POP_SIZE][G_LENGTH];
	double fitness[POP_SIZE];
	int t;
	FILE *fp;

	/* 適応度の変化を記録するファイルのオープン */
	if ((fp = fopen("result.dat", "w")) == NULL) {
		printf("Cannot open \"result.dat\"\n");
		exit(1);
	}

	/* シミュレーション条件の表示 */
	printf("個体数     : %d\n", POP_SIZE);
	printf("遺伝子長   : %d bit\n", G_LENGTH);
	printf("突然変異率 : %lf\n", M_RATE);


	init_gene(gene);              /* 遺伝子の初期化 */
	calc_fitness(gene, fitness);   /* 適応度の計算 */
	show_gene(0, gene, fitness, fp); /* 表示 */

	for (t = 1; t <= MAX_GEN; t++) {
		printf("<< 世代数 : %d >>\n", t);
		elite(gene, fitness);           /* エリート保存 */
		reproduction(gene, fitness);    /* ルーレット選択 */
		crossover(gene);               /* 単純交叉 */
		mutation(gene);                /* 突然変異 */
		calc_fitness(gene, fitness);    /* 適応度の計算 */
		show_gene(t, gene, fitness, fp);  /* 表示 */
	}

	getchar();
	fclose(fp);

	return 0;
}







