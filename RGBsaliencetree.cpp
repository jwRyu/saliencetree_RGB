#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <assert.h>
#include <time.h>
#include <opencv2/opencv.hpp>
#include <filesystem>
#include <iostream>
#include <fstream>

using namespace std;

#define OUTPUT_FNAME "C:/Users/jwryu/RUG/2018/AlphaTree/SalienceTree_colour_10rep.dat"

#define INPUTIMAGE_DIR	"C:/Users/jwryu/Google Drive/RUG/2018/AlphaTree/imgdata/Colour"
#define REPEAT 10

#define max(a,b) (a)>(b)?(a):(b)
#define min(a,b) (a)>(b)?(b):(a)
#define BOTTOM (-1)
#define false 0
#define true  1
#define MIN(a,b)  ((a<=b) ? (a) : (b))
#define MAX(a,b)  ((a>=b) ? (a) : (b))

#define CONNECTIVITY  4

typedef unsigned char uint8;
double RGBweight[3] = {0.50, 0.5, 0.5};

double MainEdgeWeight = 1.0;
double OrthogonalEdgeWeight = 1.0;

int lambda;
double omegafactor = 200000;

typedef uint8 pixel[3];
typedef unsigned long uint32;

pixel *gval=NULL, *out=NULL;

double nrmsd = 0;
uint64 memuse, max_memuse; // monitor memory usage 

inline void* Malloc(size_t size)
{
	void* pNew = malloc(size + sizeof(size_t));

	memuse += size;
	max_memuse = max(memuse, max_memuse);

	*((size_t*)pNew) = size;
	return (void*)((size_t*)pNew + 1);
}

inline void* Realloc(void* ptr, size_t size)
{
	void* pOld = (void*)((size_t*)ptr - 1);
	size_t oldsize = *((size_t*)pOld);
	void* pNew = realloc(pOld, size + sizeof(size_t));

	if (pOld != pNew)
		max_memuse = max(memuse + size, max_memuse);
	else
		max_memuse = max(memuse + size - oldsize, max_memuse);
	memuse += size - oldsize;

	*((size_t*)pNew) = size;
	return (void*)((size_t*)pNew + 1);
}

inline void Free(void* ptr)
{
	size_t size = *((size_t*)ptr - 1);
	memuse -= size;
	free((void*)((size_t*)ptr - 1));
}

typedef struct Edge{
  int p, q;
//  double alpha;
} Edge;

typedef struct{
  int maxsize;
  Edge *queue;
  uint8 *dimg;
  uint32 *bottom, *cur;
  uint32 minalpha, maxalpha;
} EdgeQueue;

EdgeQueue *EdgeQueueCreate(long maxsize, uint8 *dimg, uint32* dhist){
  EdgeQueue *newQueue = (EdgeQueue *) Malloc(sizeof(EdgeQueue));
  uint32 sum, i;
  newQueue->queue = (Edge *)Malloc(maxsize*sizeof(Edge));
  newQueue->dimg = dimg;
  newQueue->maxsize = maxsize;
  newQueue->bottom = (uint32*)Malloc(257 * sizeof(uint32));
  newQueue->cur = (uint32*)Malloc(257 * sizeof(uint32));
  newQueue->minalpha = 256;
  newQueue->maxalpha = 255;
  sum = 0;
  for (i = 0; i < 256; i++)
  {
	  newQueue->bottom[i] = newQueue->cur[i] = sum;
	  sum += dhist[i];
  }
  newQueue->cur[256] = 1;
  newQueue->bottom[256] = 0;

  return newQueue;
}


//#define EdgeQueueFront(queue)       (queue->queue + 1)
#define IsEmpty(queue)        (queue->minalpha > queue->maxalpha)

void EdgeQueueDelete(EdgeQueue *oldqueue){
  Free(oldqueue->queue);
  Free(oldqueue->bottom);
  Free(oldqueue->cur);
  Free(oldqueue->dimg);
  Free(oldqueue);
}

Edge* EdgeQueueFront(EdgeQueue *queue)
{
	return queue->queue + queue->cur[queue->minalpha] - 1;
}

void EdgeQueuePop(EdgeQueue *queue ){
  queue->cur[queue->minalpha]--;

  while (queue->cur[queue->minalpha] == queue->bottom[queue->minalpha])
	  queue->minalpha++;
}

void EdgeQueuePush(EdgeQueue *queue, int p, int q, uint8 alpha){
	uint32 idx = queue->cur[alpha]++;

	queue->queue[idx].p = p;
	queue->queue[idx].q = q;
	queue->minalpha = min(queue->minalpha, alpha);
}

typedef struct SalienceNode  
{ 
  int parent;
  int area;
  bool filtered; /* indicates whether or not the filtered value is OK */
  pixel outval;  /* output value after filtering */
  uint8 alpha;  /* alpha of flat zone */
  double sumPix[3];
  pixel minPix;
  pixel maxPix;
} SalienceNode;



typedef struct SalienceTree {
  int maxSize;
  int curSize;
  SalienceNode *node;
} SalienceTree;


SalienceTree *CreateSalienceTree(int imgsize){
  SalienceTree *tree = (SalienceTree*)Malloc(sizeof(SalienceTree));
  tree->maxSize = 2*imgsize;  /* potentially twice the number of nodes as pixels exist*/
  tree->curSize = imgsize;    /* first imgsize taken up by pixels */
  tree->node = (SalienceNode*)Malloc((tree->maxSize)*sizeof(SalienceNode));
  return tree;
}

void DeleteTree(SalienceTree *tree){
  Free(tree->node);
  Free(tree);
}


int NewSalienceNode(SalienceTree *tree, int *root, double alpha){
  SalienceNode *node =tree->node + tree->curSize;
  int result = tree->curSize;
  tree->curSize++;
  node->area = 0;
  node->alpha=alpha;
  node->parent=BOTTOM;
  root[result]=BOTTOM;
  return result;
}

void MakeSet(SalienceTree *tree, int *root, pixel *gval, int p){
  int i;
  tree->node[p].parent=BOTTOM;
  root[p]=BOTTOM;
  tree->node[p].alpha = 0.0;
  tree->node[p].area = 1;
  for (i=0;i<3;i++){
    tree->node[p].sumPix[i] = gval[p][i];
    tree->node[p].minPix[i] = gval[p][i];
    tree->node[p].maxPix[i] = gval[p][i];
  }
}

int FindRoot(int *root, int p){
  int r=p, i,j;

  while (root[r]!=BOTTOM){
    r=root[r];
  }
  i=p;
  while (i!=r){
    j=root[i];
    root[i]=r;
    i=j;
  }
  return r;
}

int FindRoot1(SalienceTree *tree, int *root, int p){
  int r=p, i,j;

  while (root[r]!=BOTTOM){
    r=root[r];
  }
  i=p;
  while (i!=r){
    j=root[i];
    root[i]=r;
    tree->node[i].parent = r;
    i=j;
  }
  return r;
}


double simpleSalience(pixel p, pixel q){
  double result=0; 
  int i;

  for (i=0;i<3;i++)
    result+= ((double)p[i] -(double)q[i])*((double)p[i] -(double)q[i]);
  return sqrt(result);
}

double WeightedSalience(pixel p, pixel q){
  double result=0; 
  int i;

  for (i=0;i<3;i++)
    result+= RGBweight[i]*((double)p[i] -(double)q[i])*((double)p[i] -(double)q[i]);
  return sqrt(result);
}


uint8 LinfNormX(pixel *img,
	int width,
	int height,
	int x, int y) {
	
	int p = width * y + x - 1, q = width * y + x;

	uint8 ret = (uint8)abs((int)img[p][0] - (int)img[q][0]);

	ret = MAX(ret, (uint8)abs((int)img[p][1] - (int)img[q][1]));
	return MAX(ret, (uint8)abs((int)img[p][2] - (int)img[q][2]));
}

uint8 LinfNormY(pixel *img,
	int width,
	int height,
	int x, int y) {

	int p = width * (y - 1) + x, q = width * y + x;

	uint8 ret = (uint8)abs((int)img[p][0] - (int)img[q][0]);

	ret = MAX(ret, (uint8)abs((int)img[p][1] - (int)img[q][1]));
	return MAX(ret, (uint8)abs((int)img[p][2] - (int)img[q][2]));
}


double EdgeStrengthX(pixel *img, 
		     int width, 
		     int height,
		     int x, int y){
  int yminus1 = y - (y > 0);
  int yplus1 = y + (y < height-1);

  double ygrad=  MIN(WeightedSalience(img[width*yminus1 + x - 1],
				      img[width*yplus1 + x - 1]),
		     WeightedSalience(img[width*yminus1 + x],
				      img[width*yplus1 + x])
		     );
  return OrthogonalEdgeWeight*ygrad + MainEdgeWeight*
    WeightedSalience(img[width*y + x - 1],
		     img[width*y + x]);
}

double EdgeStrengthY(pixel *img, 
		     int width, 
		     int height,
		     int x, int y){
  int xminus1 = x - (x > 0);
  int xplus1 = x + (x < width-1);

  double xgrad=  MIN(WeightedSalience(img[width*y + xplus1],
				      img[width*y + xminus1]),
		     WeightedSalience(img[width*(y-1) + xplus1],
				      img[width*(y-1) + xminus1])
		     );
  return OrthogonalEdgeWeight*xgrad + MainEdgeWeight*
    WeightedSalience(img[width*(y-1) + x],
		     img[width*y + x]);
}


bool IsLevelRoot(SalienceTree *tree, int i){
  int parent = tree->node[i].parent;

  if (parent==BOTTOM)
    return true;
  return (tree->node[i].alpha != tree->node[parent].alpha);
}

int LevelRoot(SalienceTree *tree, int p){
  int r=p, i,j;

  while (!IsLevelRoot(tree,r)){
    r=tree->node[r].parent;
  }
  i=p;

  while (i!=r){
    j=tree->node[i].parent;
   tree->node[i].parent=r;
    i=j;
  }
  return r;
}

# define Par(tree,p) LevelRoot(tree,tree->node[p].parent)
 

void Union(SalienceTree *tree,int *root, int p, int q){ /* p is always current pixel */
  int i;
  
  q=FindRoot1(tree,root,q);
  
  if (q!=p){
    tree->node[q].parent = p;
    root[q] = p;
    tree->node[p].area +=  tree->node[q].area;
    for (i=0;i<3;i++){
      tree->node[p].sumPix[i] += tree->node[q].sumPix[i];
      tree->node[p].minPix[i] = MIN(tree->node[p].minPix[i],tree->node[q].minPix[i]);
      tree->node[p].maxPix[i] = MAX(tree->node[p].maxPix[i],tree->node[q].maxPix[i]);
    }
  }
}


void Union2(SalienceTree *tree, int *root, int p, int q) {
	int i;
	tree->node[q].parent = p;
	root[q] = p;
	tree->node[p].area += tree->node[q].area;
	for (i = 0; i < 3; i++) {
		tree->node[p].sumPix[i] += tree->node[q].sumPix[i];
		tree->node[p].minPix[i] = MIN(tree->node[p].minPix[i], tree->node[q].minPix[i]);
		tree->node[p].maxPix[i] = MAX(tree->node[p].maxPix[i], tree->node[q].maxPix[i]);
	}
}

void Union3(SalienceTree *tree, int *root, int p, int q, int r) {
	int i;
	tree->node[p].parent = r;
	tree->node[q].parent = r;
	root[p] = r;
	root[q] = r;
	tree->node[r].area = tree->node[p].area + tree->node[q].area;
	for (i = 0; i < 3; i++) {
		tree->node[r].sumPix[i] = tree->node[p].sumPix[i] + tree->node[q].sumPix[i];
		tree->node[r].minPix[i] = MIN(tree->node[p].minPix[i], tree->node[q].minPix[i]);
		tree->node[r].maxPix[i] = MAX(tree->node[p].maxPix[i], tree->node[q].maxPix[i]);
	}
}

void compute_dhist(uint32 *dhist, uint8 *dimg, pixel *img, int width, int height) {
	/* pre: tree has been created with imgsize= width*height
	queue initialized accordingly;
	*/
	uint32 imgsize = width * height;
	uint32 p, x, y;
	uint32 dimgidx;
	uint8 edgeSalience;

	dimgidx = 3;
	for (x = 1; x < width; x++) {
		edgeSalience = LinfNormX(img, width, height, x, 0);
		dimg[dimgidx] = edgeSalience;
		dhist[edgeSalience]++;
		dimgidx += 2;
	}
	dimgidx--;

	for (y = 1; y < height; y++) {
		p = y * width;
		edgeSalience = LinfNormY(img, width, height, 0, y);
		dimg[dimgidx] = edgeSalience;
		dhist[edgeSalience]++;
		dimgidx += 2;

		p++;
		for (x = 1; x < width; x++, p++) {
			edgeSalience = LinfNormY(img, width, height, x, y);
			dimg[dimgidx++] = edgeSalience;
			dhist[edgeSalience]++;


			edgeSalience = LinfNormX(img, width, height, x, y);
			dimg[dimgidx++] = edgeSalience;
			dhist[edgeSalience]++;
		}
	}
}

void Phase1(SalienceTree *tree, EdgeQueue *queue, int *root,
	pixel *img,	int width, int height, double lambdamin) {
	/* pre: tree has been created with imgsize= width*height
	queue initialized accordingly;
	*/
	uint32 imgsize = width * height;
	uint32 p, x, y;
	uint32 dimgidx;
	uint8 *dimg = queue->dimg;

	uint8 edgeSalience;

	MakeSet(tree, root, img, 0);

	dimgidx = 3;
	for (x=1;x<width;x++){
		MakeSet(tree, root, img, x);
		edgeSalience = dimg[dimgidx];
		dimgidx += 2;
		if (edgeSalience == 0)
			Union(tree, root, x, x-1);
		else
			EdgeQueuePush(queue,x,x-1,edgeSalience);
	}
	dimgidx--;

	for (y = 1; y < height; y++){
		p = y * width;
		MakeSet(tree, root, img, p);
		edgeSalience = dimg[dimgidx];
		dimgidx += 2;

		if (edgeSalience == 0)
			Union(tree, root, p, p-width);
		else
			EdgeQueuePush(queue,p, p-width,edgeSalience);

		p++;  
		for (x=1;x<width;x++,p++){
			MakeSet(tree, root, img, p);
			edgeSalience = dimg[dimgidx++];

			if (edgeSalience == 0)
				Union(tree, root, p, p-width);
			else
				EdgeQueuePush(queue,p, p-width,edgeSalience);
		
			edgeSalience = dimg[dimgidx++];

			if (edgeSalience == 0)
				Union(tree, root, p, p-1);
			else
				EdgeQueuePush(queue,p, p-1,edgeSalience);
		}    

	}

}

void GetAncestors(SalienceTree *tree, int *root, int *p, int*q){
	int temp;
	*p = LevelRoot(tree, *p);
	*q = LevelRoot(tree, *q);
	if (*p<*q){
		temp = *p;
		*p=*q;
		*q=temp;
	}
	while ((*p!=*q) && (root[*p]!=BOTTOM) &&(root[*q]!=BOTTOM)){
		*q=root[*q];
		if (*p<*q){
			temp = *p;
			*p=*q;
			*q=temp;
		}
	}
	if (root[*p]==BOTTOM) {
		*q=FindRoot(root,*q);
	} else if (root[*q]==BOTTOM){
		*p=FindRoot(root,*p);
	}
}

void Phase2(SalienceTree *tree, EdgeQueue *queue, int *root, 
	pixel *img, 
	int width, int height){
	Edge *currentEdge;
	int v1, v2, temp, r;
	uint8 oldalpha=0, alpha12; 
	while (!IsEmpty(queue)){
		currentEdge = EdgeQueueFront(queue);
		v1 =currentEdge->p;
		v2 =currentEdge->q;
		GetAncestors(tree,root,&v1,&v2);
		alpha12 = queue->minalpha;

		EdgeQueuePop(queue);
		if (v1!=v2) {
			if (v1 < v2) {
				temp=v1;
				v1=v2;
				v2=temp;
			}
			if (tree->node[v1].alpha < alpha12){
				r =NewSalienceNode(tree,root, alpha12);
				Union3(tree,root,v1,v2,r);
			} else {
				Union2(tree,root,v1,v2);
			}
		}
		oldalpha=alpha12;
	}
}

SalienceTree *MakeSalienceTree(pixel *img,  
	int width, int height, int channel, double lambdamin){
	int imgsize = width*height;
	EdgeQueue *queue;
	int *root = (int*)Malloc(imgsize*2*sizeof(int));
	SalienceTree *tree;
	uint32 *dhist;
	uint8 *dimg;

	dhist = (uint32*)Malloc(256 * sizeof(uint32));
	dimg = (uint8*)Malloc(imgsize * 2 * sizeof(uint8));
	memset(dhist, 0, 256 * sizeof(uint32));

	compute_dhist(dhist, dimg, img, width, height);
	queue = EdgeQueueCreate((CONNECTIVITY / 2)*imgsize, dimg, dhist);
	tree = CreateSalienceTree(imgsize);
	assert(tree!=NULL);
	assert(tree->node!=NULL);
	//fprintf(stderr,"Phase1 started\n");
	Phase1(tree,queue,root,img,width,height,lambdamin);
	//fprintf(stderr,"Phase2 started\n");
	Phase2(tree,queue,root,img,width,height);
	//fprintf(stderr,"Phase2 done\n");
	EdgeQueueDelete(queue);
	Free(root);

	Free(dhist);
	return tree;
}

void SalienceTreeAreaFilter(SalienceTree *tree, pixel *out, int lambda){
	int i,j, imgsize = tree->maxSize/2;
	if (lambda <= imgsize){
		for (j=0;j<3;j++){
			tree->node[tree->curSize-1].outval[j] =
				tree->node[tree->curSize-1].sumPix[j]/tree->node[tree->curSize-1].area;
		}
		for (i = tree->curSize-2; i>=0; i--){

			if (IsLevelRoot(tree,i) && (tree->node[i].area>=lambda)){
				for (j=0;j<3;j++)
					tree->node[i].outval[j] = tree->node[i].sumPix[j]/tree->node[i].area;
			} else {
				for (j=0;j<3;j++)
					tree->node[i].outval[j] = tree->node[tree->node[i].parent].outval[j];
			}
		}
	} else {
		for (i = tree->curSize-1; i>=0; i--){
			for (j=0;j<3;j++)
				tree->node[i].outval[j] = 0;
		}
	}
	for (i = 0; i < imgsize; i++)
		for (j=0;j<3;j++)
			out[i][j] = tree->node[i].outval[j];
}

#define NodeSalience(tree, p) (tree->node[Par(tree,p)].alpha)

void SalienceTreeSalienceFilter(SalienceTree *tree, pixel *out, double lambda){
	int i,j, imgsize = tree->maxSize/2;
	if (lambda <= tree->node[tree->curSize-1].alpha){
		for (j=0;j<3;j++){
			tree->node[tree->curSize-1].outval[j] =
				tree->node[tree->curSize-1].sumPix[j]/tree->node[tree->curSize-1].area;
		}
		for (i = tree->curSize-2; i>=0; i--){

			if (IsLevelRoot(tree,i) && (NodeSalience(tree,i)>=lambda)){
				for (j=0;j<3;j++)
					tree->node[i].outval[j] = tree->node[i].sumPix[j]/tree->node[i].area;
			} else {
				for (j=0;j<3;j++)
					tree->node[i].outval[j] = tree->node[tree->node[i].parent].outval[j];
			}
		}
	} else {
		for (i = tree->curSize-1; i>=0; i--){
			for (j=0;j<3;j++)
				tree->node[i].outval[j] = 0;
		}
	}
	for (i = 0; i < imgsize; i++)
		for (j=0;j<3;j++)
			out[i][j] = tree->node[i].outval[j];

}


/*

short ImagePPMAsciiRead(char *fname)
{
   FILE *infile;
   ulong i,j;
   int c;

   infile = fopen(fname, "r");
   if (infile==NULL) {
      fprintf (stderr, "Error: Can't read the ASCII file: %s !", fname);
      return(0);
   }
   fscanf(infile, "P3\n");
   while ((c=fgetc(infile)) == '#')
      while ((c=fgetc(infile)) != '\n');
   ungetc(c, infile);
   fscanf(infile, "%d %d\n255\n", &width, &height);
   size = width*height;

   gval = malloc(size*sizeof(Pixel));
   if (gval==NULL) {
      fprintf (stderr, "Out of memory!");
      fclose(infile);
      return(0);
   }
   for (i=0; i<size; i++)
   {
     for (j=0;j<3;j++){
       fscanf(infile, "%d", &c);
       gval[i][j] = c;
     }
   }
   fclose(infile);
   return(1);
}


short ImagePPMBinRead(char *fname)
{
   FILE *infile;
   int c, i;

   infile = fopen(fname, "rb");
   if (infile==NULL) {
      fprintf (stderr, "Error: Can't read the binary file: %s !", fname);
      return(0);
   }
   fscanf(infile, "P6\n");
   while ((c=fgetc(infile)) == '#')
      while ((c=fgetc(infile)) != '\n');
   ungetc(c, infile);
   fscanf(infile, "%d %d\n255\n", &width, &height);
   size = width*height;

   gval = malloc(size*sizeof(Pixel));
   if (gval==NULL) {
     fprintf (stderr, "Out of memory!");
     fclose(infile);
     return(0);
   } 
   fread(gval, sizeof(Pixel), size, infile);

   fclose(infile);
   return(1);
} 


short ImagePPMRead(char *fname)
{
   FILE *infile;
   char id[4];

   infile = fopen(fname, "r");
   if (infile==NULL) {
      fprintf (stderr, "Error: Can't read the image: %s !", fname);
      return(0);
   }
   fscanf(infile, "%3s", id);
   fclose(infile);
   if (strcmp(id, "P3")==0) return(ImagePPMAsciiRead(fname));
   else if (strcmp(id, "P6")==0) return(ImagePPMBinRead(fname));
   else {
     fprintf (stderr, "Unknown type of the image!");
     return(0);
   }
} 

int ImagePPMBinWrite(char *fname)
{
   FILE *outfile;

   outfile = fopen(fname, "wb");
   if (outfile==NULL) {
      fprintf (stderr, "Error: Can't write the image: %s !", fname);
      return(-1);
   }
   fprintf(outfile, "P6\n%d %d\n255\n", width, height);

   fwrite(out,sizeof(Pixel) , (size_t)(size), outfile);

   fclose(outfile);
   return(0);
} 
*/

// reshape 3d image matrix (ch,x,y) -> (x,y,ch)
void imreshape(uint8* dst, uint8* src, uint32 height, uint32 width)
{
	uint32 ch, i, dstidx, srcidx;
	srcidx = 0;
	for (ch = 0; ch < 3; ch++)
	{
		dstidx = ch;
		for (i = 0; i < height * width; i++)
		{
			dst[dstidx] = src[srcidx++];
			dstidx += 3;
		}
	}
}

int main(int argc, char *argv[]) {
	uint32 width, height, channel;
	double lambda;
	SalienceTree *tree;
	uint32 cnt = 0;
	pixel *img; 	
	ofstream f;
	ifstream fcheck;
	uint32 contidx;
	char in;
	contidx = 0;
	//	f.open("C:/Users/jwryu/RUG/2018/AlphaTree/AlphaTree_grey_Exp.dat", std::ofstream::app);
	fcheck.open(OUTPUT_FNAME);
	if (fcheck.good())
	{
		cout << "Output file \"" << OUTPUT_FNAME << "\" already exists. Overwrite? (y/n/a)";
		cin >> in;
		if (in == 'a')
		{
			f.open(OUTPUT_FNAME, std::ofstream::app);
			cout << "Start from : ";
			cin >> contidx;
		}
		else if (in == 'y')
			f.open(OUTPUT_FNAME);
		else
			exit(-1);
	}
	else
		f.open(OUTPUT_FNAME);

	lambda = 3;
	cnt = 0;
	std::string path = INPUTIMAGE_DIR;
	for (auto & p : std::experimental::filesystem::directory_iterator(path))
	{
		if (++cnt < contidx)
		{
			cout << cnt << ": " << p << ' ' << endl;
			continue;
		}
		cv::String str1(p.path().string().c_str());
		cv::Mat cvimg = imread(str1, cv::IMREAD_ANYCOLOR);

		height = cvimg.rows;
		width = cvimg.cols;
		channel = cvimg.channels();

		cout << cnt << ": " << str1 << ' ' << height << 'x' << width << endl;

		if (channel != 3)
		{
			cout << "input should be a 3-ch image" << endl;
			getc(stdin);
			exit(-1);
		}

		img = (pixel*)malloc(height * width * sizeof(pixel));
		imreshape((uint8*)img, cvimg.data, height, width);
		cvimg.release();
		str1.clear();
		
		double runtime, minruntime;
		for (int testrep = 0; testrep < REPEAT; testrep++)
		{
			memuse = max_memuse = 0;
			auto wcts = std::chrono::system_clock::now();

			tree = MakeSalienceTree((pixel*)img, width, height, channel, (double)lambda);

			std::chrono::duration<double> wctduration = (std::chrono::system_clock::now() - wcts);
			runtime = wctduration.count();
			minruntime = testrep == 0 ? runtime : min(runtime, minruntime);

			if (testrep < (REPEAT - 1))
				DeleteTree(tree);
		}
		f << p.path().string().c_str() << '\t' << height << '\t' << width << '\t' << max_memuse << '\t' << nrmsd << '\t' << tree->maxSize << '\t' << tree->curSize << '\t' << minruntime << endl;

		cout << "Time Elapsed: " << minruntime << endl;

		free(img);
		DeleteTree(tree);
	}
   return(0);
} /* main */
