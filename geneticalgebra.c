/********************************************************************
��`�I�A���S���Y�� 
********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>

#define POP_SIZE 5   /* �̐� (�K����ɐݒ�) */
#define G_LENGTH 10  /* �̂̈�`�q�^�̃r�b�g�� */
#define MAX_GEN  20  /* ���㐔 */
#define M_RATE   0.1 /* �ˑR�ψٗ� (0?1) */ 
#define MAX_Fitt 10 /* �ō��K����*/ 

/********************************************************************
�����̔��� (Seed�̌���)
********************************************************************/
void init_rnd()
{
	srand((unsigned int)time(NULL));
}

/********************************************************************
�����̔��� (0?1�̗���)
********************************************************************/
double Random()
{
	return((double)rand() / RAND_MAX);
}

/********************************************************************
��`�q�̏�����
���� gene[p][i] : ��`�qp��i�Ԗڂ̐���
********************************************************************/
void init_gene(int gene[POP_SIZE][G_LENGTH])
{
	int p, i;

	/* �����̏����� */
	init_rnd();

	/* ��`�q��������  0?1�̗����𔭐����A0.5�ȏ�Ȃ�1
	0.5�����Ȃ�0 */
	printf("<< �����̌Q >>\n");

	for (i = 0; i < POP_SIZE; i++) {//�����̒ǉ��Ǝd��
		for (p = 0; p < G_LENGTH; p++) {
			gene[i][p] = rand()%2;
			printf("%d%d\n", gene[i][p]);
		}
	}

	for (i = 0; i < POP_SIZE; i++) {//���l�̐���
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
�K���x�̌v�Z
���� gene[p][i] : ��`�qp��i�Ԗڂ̐���
fitness[p] : ��`�qp�̓K���x
********************************************************************/
void calc_fitness(int gene[POP_SIZE][G_LENGTH], double fitness[POP_SIZE])
{
	int p, i;

	/* �K���x�̌v�Z �O����5bit��0�̐� �㔼��5bit��1�̐� */

	for (i = 0; i < POP_SIZE; i++) {
		fitness[i] = 0;//�������т̒ǉ�
		for (p = 0; p < G_LENGTH; p++) {
			if (p < 5 && gene[i][p] == 0) {//�K���x�̌v�Z
				fitness[i]++;
			}
			else if (p <= 10 && gene[i][p] == 1) {
				fitness[i]++;
			}
		}
	}
}

/**********************************************************************
��`�q�̕\�� & �ő�K���x�E���ϓK���x�̌v�Z & �t�@�C���ւ̏����o��
���� t          : ���㐔
gene[p][i] : ��`�qp��i�Ԗڂ̐���
fitness[p] : ��`�qp�̓K���x
*fp        : �t�@�C���|�C���^
**********************************************************************/
void show_gene(int t, int gene[POP_SIZE][G_LENGTH],
	double fitness[POP_SIZE],
	FILE *fp)
{
	int p, i;
	double avg_fit=0;  /* ���ϓK���x  */
	double max_fit=0;  /* �ő�K���x  */

					 /* �̂̒l�A�K���x�̕\�� */
					 /* ���ρE�ő�K���x�̌v�Z */

	printf("�K���x�\��");
	for (i = 0; i < POP_SIZE; i++) {
		for (p = 0; p < G_LENGTH; p++) {
			printf("%d\n", fitness[i]);
			if (max_fit<fitness[i]) {
				max_fit = gene[i][p];//�ő�K���̌̌Q��ۑ�
			}
			avg_fit += gene[i][p];//���Z
		}
	}
	avg_fit += avg_fit / POP_SIZE;//���ϒl�̌v�Z
	
	
					 /* ���ρE�ő�K���x�̕\�� */
	printf("���ϓK���x : %lf\n", avg_fit);
	printf("�ő�K���x : %lf\n", max_fit);

	/* ���ρE�ő�K���x���t�@�C���ɏ������� */
	fprintf(fp, "%d %lf %lf\n", t, avg_fit, max_fit);
}

/**********************************************************************
�̔ԍ� p1 �� p2 �̓K���x�ƈ�`�q������
���� p1, p2     : ��`�q�̔ԍ�
gene[p][i] : ��`�qp��i�Ԗڂ̐���
fitness[p] : ��`�qp�̓K���x
**********************************************************************/
void swap_gene(int p1, int p2, int gene[POP_SIZE][G_LENGTH],
	double fitness[POP_SIZE])
{
	int tmp;
	double f;
	int i,p;

	/* ��`�q�^�̌��� (��`�qp1�ƈ�`�qp2�̒l�����ւ���) */

	while (true) {
		//���̐���̌̌Q�̑S�Ă̌̂̐��сi�I�񂾕i���̉��l�̍��v�j���Z�o����
		for (i = 0; i < POP_SIZE; i++) {
			fitness[i] = 0;//�����̂̐���
			for (p = 0; p < G_LENGTH; p++) {
				fitness[i] += gene[i][p];
			}
			if (fitness[i] > MAX_GEN) {
				fitness[i] = 0.0;
			}
		}
	}

	/* �K���x�̌��� (��`�qp1�ƈ�`�qp2�̓K���x�̒l�����ւ���) */
}

/**********************************************************************
�̔ԍ� p1 �̓K���x�ƈ�`�q�^�� p2 �ɃR�s�[
���� p1, p2     : ��`�q�̔ԍ�
gene[p][i] : ��`�qp��i�Ԗڂ̐���
fitness[p] : ��`�qp�̓K���x
**********************************************************************/
void copy_gene(int p1, int p2, int gene[POP_SIZE][G_LENGTH],
	double fitness[POP_SIZE])
{
	int i;

	/* ��`�q�̃R�s�[ (��`�qp1����`�qp2�ɃR�s�[����) */
	/* �K���x�̃R�s�[ (��`�qp1�̓K���x����`�qp2�̓K���x�ɃR�s�[����)*/
	while (true) {
		//���̐���̌̌Q�̑S�Ă̌̂̐��сi�I�񂾕i���̉��l�̍��v�j���Z�o����
		for (i = 0; i < POP_SIZE; i++) {
			fitness[i] = 0;//�����̂̐���
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
�G���[�g�ۑ�
(�ŏ��K���x�̌̂ɍő�K���x�̌̂̃f�[�^���R�s�[)
���� gene[p][i] : ��`�qp��i�Ԗڂ̐���
fitness[p] : ��`�qp�̓K���x
**********************************************************************/
void elite(int gene[POP_SIZE][G_LENGTH], double fitness[POP_SIZE])
{
	int p, i;
	double max_fitness = fitness[0];
	double min_fitness = fitness[0];
	int max_p = 0;
	int min_p = 0;
	int M_F

	/* �ő�K���x�̌�(max_p)�ƍŏ��K���x�̌�(min_p)�������� */

	//�ő�K���x
	for (i = 0; i < POP_SIZE; i++) {
		if (max_fitness < fitness[i]) {
			max_fitness < fitness[i];//���т��\�[�g���čő�l��ۑ�
			M_F = i;//�G���[�g�I�����ꂽ�ԍ����L�^
		}
	}
	//�ŏ��K���x
	for (i = 0; i <POP_SIZE; i++) {
			if (max_fitness < fitness[i]) {//�ŗD�G��`�q�̕ۑ�
				min_p = gene[M_F][i];
			}
	}


	/* �ŏ��K���x�̌̂ɍő�K���x�̌̂��R�s�[ */
	copy_gene(max_p, min_p, gene, fitness);
	/* �ő�K���x�̌̂�0�ԖڂɈړ� */
	swap_gene(0, max_p, gene, fitness);
}

/**********************************************************************
���[���b�g�I��
���� gene[p][i] : ��`�qp��i�Ԗڂ̐���
fitness[p] : ��`�qp�̓K���x
**********************************************************************/
void reproduction(int gene[POP_SIZE][G_LENGTH], double fitness[POP_SIZE])
{
	double sum_of_fitness; /* �̂̓K���x�̑��a */
	double border;         /* ���[���b�g��̌̊Ԃ̋��E */
	double r;              /* ���[���b�g��̑I���ʒu */
	int p, i;               /* �I�΂ꂽ�̂̔ԍ� */
	int num;               /* 0 <= num <= POP_SIZE-1 */
	int new_gene[POP_SIZE][G_LENGTH];

	/* ���[���b�g��1���� sum_of_fitness �����߂� */
	sum_of_fitness = 0.0;
	for (p = 0; p<POP_SIZE; p++) {
		sum_of_fitness += fitness[p];
	}

	/* ���[���b�g�� POP_SIZE �񂾂��񂵂Ď�����̌̂�I�� */
	for (p = 1; p<POP_SIZE; p++) {
		/* ���[���b�g���񂵂ďꏊ��I��
		r : �I�΂ꂽ�ʒu (0 <= r <= sum_of_fitness) */
		r = sum_of_fitness * Random();
		/* �I�΂ꂽ�ꏊ�ɊY������̂����Ԃ����ׂ�
		num : �I�΂ꂽ�̂̔ԍ� (0 <= num <= POP_SIZE-1) */
		num = 0;
		border = fitness[0]; /* �̊Ԃ̋��E */
		while (border<r) {
			num++;
			border += fitness[num];
		}

		/* ��`�q�̑�� */
		for (i = 0; i<G_LENGTH; i++) {
			new_gene[p][i] = gene[num][i];
		}
	}

	/* ��`�q�̃R�s�[ */
	for (p = 1; p<POP_SIZE; p++) {
		for (i = 0; i<G_LENGTH; i++) {
			gene[p][i] = new_gene[p][i];
		}
	}
}

/**********************************************************************
��_����
���� gene[p][i] : ��`�qp��i�Ԗڂ̐���
**********************************************************************/
void crossover(int gene[POP_SIZE][G_LENGTH])
{
	int gene1[G_LENGTH]; /* �e1�̈�`�q�^ */
	int gene2[G_LENGTH]; /* �e2�̈�`�q�^ */
	int i, j;
	int c_pos;   /* �����ʒu (1 <= c_pos <= G_LENGTH-1) */

				 /* �����ʒu��1?G_LENGTH-1�͈̔͂Ń����_���Ɍ��߁A
				 ������������ւ���B
				 gene[1]��gene[2],  gene[3]��gene[4] ... �̂悤�ɐe�ɂ��� */

	for (i = 0; i < POP_SIZE;i+=2) {
		c_pos = rand()%G_LENGTH;//�؏��_�������_���Ɏw��
		for (j = 0; j < G_LENGTH; j++) {
			gene2[i] = gene[i][j];
			gene[i][j] = gene[i + 1][j];//��x���̒l�Ɉ�`�q�𓦂����Ă���đ�����邱�ƂŌ���������
			gene[i + 1][j] = gene2[i];
		}
	}

}

/**********************************************************************
�ˑR�ψ�
���� gene[p][i] : ��`�qp��i�Ԗڂ̐���
**********************************************************************/
void mutation(int gene[POP_SIZE][G_LENGTH])
{
	int p, i;
	double random;

	/* 0?1�̗����𔭐������A���̒l�� M_RATE �ȉ��Ȃ��
	��`�q�̒l�������_���ɕς��� (0�Ȃ��1�A1�Ȃ��0) */

	if (i = 0; i < POP_SIZE; i++) {
		for (p = 0; p < G_LENGTH; p++) {
			random = (double)rand() / 101 / 100;//�ˑR�ψٗ�
			if (random < M_RATE) {
				gene[i][p] = rand() % 2;//�����_����1��0����
			}
		}
	}

}

/**********************************************************************
���C���v���O����
**********************************************************************/
int main(int argc, char *argv[])
{
	int gene[POP_SIZE][G_LENGTH];
	double fitness[POP_SIZE];
	int t;
	FILE *fp;

	/* �K���x�̕ω����L�^����t�@�C���̃I�[�v�� */
	if ((fp = fopen("result.dat", "w")) == NULL) {
		printf("Cannot open \"result.dat\"\n");
		exit(1);
	}

	/* �V�~�����[�V���������̕\�� */
	printf("�̐�     : %d\n", POP_SIZE);
	printf("��`�q��   : %d bit\n", G_LENGTH);
	printf("�ˑR�ψٗ� : %lf\n", M_RATE);


	init_gene(gene);              /* ��`�q�̏����� */
	calc_fitness(gene, fitness);   /* �K���x�̌v�Z */
	show_gene(0, gene, fitness, fp); /* �\�� */

	for (t = 1; t <= MAX_GEN; t++) {
		printf("<< ���㐔 : %d >>\n", t);
		elite(gene, fitness);           /* �G���[�g�ۑ� */
		reproduction(gene, fitness);    /* ���[���b�g�I�� */
		crossover(gene);               /* �P������ */
		mutation(gene);                /* �ˑR�ψ� */
		calc_fitness(gene, fitness);    /* �K���x�̌v�Z */
		show_gene(t, gene, fitness, fp);  /* �\�� */
	}

	getchar();
	fclose(fp);

	return 0;
}







