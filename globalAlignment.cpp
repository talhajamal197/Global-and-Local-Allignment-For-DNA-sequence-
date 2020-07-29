#include<iostream>
#include<string>
using namespace std;
struct Node
{   int x;
    int y;
	Node *next;
	Node *prev;
	char ch1;
	char ch2;
	int which;
};
void loadProgram()
{
	cout<<"\n\n\n\n\n\t\t\t\t\t";
	cout<<"LOADING!\n\t\t\t\t";
	char a=219;
	for(int i=0;i<25;i++)
	{
		for(int j=0;j<10e5;j++)
		{
		}
		cout<<a;
	}
	system("cls");
}
string seq1="ATCGTCGAATCGTCGAATCGTCGAA";//taken on x axis
string seq2="TCGGGTACATTCGGGTACATT";

int rows,cols;
int **scoringMatrix;
struct Node **ptr;
int gap=-2;
int match=1;
int misMatch=-1;

int findMax(int x,int y,int *max)
{   int which=1;
	int up=scoringMatrix[x][y-1],diagonal=scoringMatrix[x-1][y-1],left=scoringMatrix[x-1][y];
	*max=up+gap;
	if(left+gap>=*max)
	{
		*max=left+gap;
		which=3;
	}
	//check if match or mismatch
	char a=seq1[y-1];
	char b=seq2[x-1];
	if(a==b&&diagonal+match>=*max)
	{
		*max=diagonal+match;
		which=2;
	}
	else if(diagonal+misMatch>=*max)
	{
	    *max=diagonal+misMatch;
		which=2;
	}
	
	return which;
}
void backTrack(struct Node *curr)
{   
    struct Node *temp=curr;
    struct Node *prev=0;
while(curr->prev!=0)
{   prev=curr;
	curr->prev->next=curr;
	curr=curr->prev;
}
temp=prev;
while(prev!=0)
{
	cout<<prev->ch1<<' ';
	prev=prev->next;
}
prev=temp;
cout<<endl;
while(prev!=0)
{
	cout<<prev->ch2<<' ';
	prev=prev->next;
}
}

void fillMatrix()
{   int max=0;
    int which=0; 
    //column wise filling
	for(int i=1;i<rows;i++)
	{
		for(int j=1;j<cols;j++)
		{
		  	which=findMax(i,j,&max);
			  scoringMatrix[i][j]=max;
			  if(which==1)
			  { //returns 1 if max is from up,2 for diagonal,3 for left
			  	ptr[i][j].prev=&ptr[i][j-1];
			  	ptr[i][j].ch1=seq1[j-1];//seq1
			  ptr[i][j].ch2='-';//seq2
			  char t1=ptr[i][j-1].ch1;//prev char 1 
			  char t2=ptr[i][j-1].ch2;//prev char 2
			  
			  
			  }
			  if(which==2)
			  {
			  	ptr[i][j].prev=&ptr[i-1][j-1];
			    ptr[i][j].ch1=seq1[j-1];//seq1
			    ptr[i][j].ch2=seq2[i-1];//seq2
			    char t1=ptr[i-1][j-1].ch1;
			  char t2=ptr[i-1][j-1].ch2;
			  }
			  if(which==3)
			  {
			ptr[i][j].prev=&ptr[i-1][j];
			  ptr[i][j].ch1='-';//seq1
			  ptr[i][j].ch2=seq2[j];//seq2
			  char t1=ptr[i-1][j].ch1;
			  char t2=ptr[i-1][j].ch2;
			  }
			  ptr[i][j].x=i;
			  ptr[i][j].y=j;
			  //string seq1="TTCA";
              //string seq2="ACT";
              
			  ptr[i][j].which=which;
			  
		}
		
	}
	
	backTrack(&ptr[rows-1][cols-1]);
	cout<<endl;
}
int main()
{   
    
    
    rows=seq2.length()+1;
    cols=seq1.length()+1;
    //Allocate space for row pointers +1 for gap
	scoringMatrix=new int*[rows];
	//creating array of nodes
	ptr=new struct Node*[rows];
	//Allocate space for col pointers +1 for gap
	for(int i=0;i<rows;i++)
	{
		scoringMatrix[i]=new int[cols];
		ptr[i]=new struct Node[cols];
	}
	//Initializing the Matrix with 0
	for(int i=0;i<rows;i++)
	{
		for(int j=0;j<cols;j++)
		{
			scoringMatrix[i][j]=0;
		}
	}
	
	
	ptr[0][0].ch1='T';
	ptr[0][0].ch2='C';
	//Fill 1st row and col with -2 incrementaion
	for(int i=0;i<1;i++)
	{   //row wise or first col.
		for(int j=1;j<cols;j++)
		{   ptr[i][j].prev=&ptr[i][j-1];
		    char a=ptr[i][j-1].ch1;
		    char b=ptr[i][j-1].ch2;
		    ptr[i][j].x=i;
		    ptr[i][j].y=j;
		    ptr[i][j].which=1;
		    ptr[i][j].ch2='-';//-
			ptr[i][j].ch1=seq1[j-1];//seq2
			scoringMatrix[i][j]+=scoringMatrix[i][j-1]-2;
		}
	}
	
	//col wise or first row.
	for(int i=0;i<1;i++)
	{
		for(int j=1;j<rows;j++)
		{   ptr[j][i].prev=&ptr[j-1][i];
		    char a=ptr[j-1][i].ch1;
		    char b=ptr[j-1][i].ch2;
		    ptr[j][i].x=j;
		    ptr[j][i].y=i;
		    ptr[j][i].which=3;
		    ptr[j][i].ch2=seq2[j-1];//-
			ptr[j][i].ch1='-';//seq2
			scoringMatrix[j][i]+=scoringMatrix[j-1][i]-2;
		}
	}
	
	int max=0;
	//scoringMatrix
	loadProgram();
//	cout<<"\nEnter DNA Sequence 1:";
//	cin>>seq1;
//	cout<<"\nEnter DNA Sequence 2:";
//	cin>>seq2;
	fillMatrix();
	
    //Display the matrix
    cout<<"Scoring Matrix"<<endl;
	for(int i=0;i<rows;i++)
	{
		for(int j=0;j<cols;j++)
		{
			cout<<scoringMatrix[i][j]<<' ';
		}
		cout<<'\n';
	}

//for(int i=0;i<rows;i++)
//{
//	for(int j=0;j<cols;j++)
//	{   cout<<"Current-Address:"<<&ptr[i][j]<<endl;
//	    cout<<"Previous-Address:"<<ptr[i][j].prev<<endl;
//		cout<<(ptr[i][j].x)<<','<<ptr[i][j].y<<' '<<"which:"<<ptr[i][j].which<<endl;
//			
//	}
//		
//}
}

