#ifndef MEMORY_H
#define MEMORY_H

#include<iostream>
#include<iomanip>

template <class type>
class matrix4D{
	public:
		int d1,d2,d3,d4;
		type**** memory;

		matrix4D();
		matrix4D(int,int,int,int);
		matrix4D(int,int,int,int,type);
		~matrix4D();

		/*********************
			Retrieve memory
		*********************/
		inline type value(int,int,int,int);
		inline void print(int,int,int,int);

		/*********************
			Alter memory
		*********************/
		void resize(int,int,int,int,bool);
		inline void set(int,int,int,int,type);
};

template <class type>
class storage{
	public:
		int rows, *cols;
		type** memory;

		storage();
		storage(int,int*);
		~storage();

		/*********************
			Retrieve memory
		*********************/
		inline type value(int,int);
		inline void print(int,int);

		/*********************
			Alter memory
		*********************/
		void resize(int,int*,bool);
		inline void set(int,int,type);
};

template <class type>
class matrix{
	public:
		int rows, cols;
		type** memory;

		matrix();
		matrix(int,int);
		matrix(int,int,type);
		~matrix();

		/*********************
			Retrieve memory
		*********************/
		inline type value(int,int);
		inline void print(int,int);

		/*********************
			Alter memory
		*********************/
		void resize(int,int,bool);
		inline void set(int,int,type);
};

template <class type>
class array{
	public:
		int length;
		type* memory;
		array();
		array(int);
		array(int,type);
		~array();

		/*********************
			Retrieve memory
		*********************/
		inline type value(int);
		inline void print(int);

		/*********************
			Alter memory
		*********************/
		void resize(int,bool);
		inline void set(int,type);
};


/*********************
	con/destuctors
*********************/
template <class type>
matrix4D<type>::matrix4D(){
	d1=0; d2=0; d3=0; d4=0;
	memory = nullptr;
}
template <class type>
matrix4D<type>::matrix4D(int x, int y, int z, int d){
	int i,j,k;
	d1=x; d2=y; d3=z; d4=d;
	if(d1==0||d2==0||d3==0||d4==0){
		memory = nullptr;
	}
	else{
		memory = new type***[d1];
		for(i=0;i<d1;i++){
			memory[i] = new type**[d2];
			for(j=0;j<d2;j++){
				memory[i][j] = new type*[d3];
				for(k=0;k<d3;k++){
					memory[i][j][k] = new type[d4];
				}
			}
		}
	}
}
template <class type>
matrix4D<type>::matrix4D(int x, int y, int z, int d, type a){
	int i,j,k,l;
	d1=x; d2=y; d3=z; d4=d;
	if(d1==0||d2==0||d3==0||d4==0){
		memory = nullptr;
	}
	else{
		memory = new type***[d1];
		for(i=0;i<d1;i++){
			memory[i] = new type**[d2];
			for(j=0;j<d2;j++){
				memory[i][j] = new type*[d3];
				for(k=0;k<d3;k++){
					memory[i][j][k] = new type[d4];
				}
			}
		}
		for(i=0;i<d1;i++){
			for(j=0;j<d2;j++){
				for(k=0;k<d3;k++){
					for(l=0;l<d4;l++){
						memory[i][j][k][l] = a;
					}
				}
			}
		}
	}
}
template <class type>
matrix4D<type>::~matrix4D(){
	int i,j,k;
	for(i=0;i<d1;i++){
		for(j=0;j<d2;j++){
			for(k=0;k<d3;k++){
				delete[] memory[i][j][k];
			}
			delete[] memory[i][j];
		}
		delete[] memory[i];
	}
	delete[] memory;
	d1=0; d2=0; d3=0; d4=0;
	memory = nullptr;
}

template <class type>
storage<type>::storage(){
	rows = 0; cols = nullptr;
	memory = nullptr;
}
template <class type>
storage<type>::storage(int x, int*y){
	int i;
	rows = x;
	cols = new int[rows];
	for(i=0;i<rows;i++){
		cols[i] = y[i];
	}
	if(rows==0||cols==nullptr){
		memory = nullptr;
	}
	else{
		memory = new type*[rows];
		for(i=0;i<rows;i++){
			memory[i] = new type[cols[i]];
		}
	}
}
template <class type>
storage<type>::~storage(){
	int i;
	for(i=0;i<rows;i++){
		delete[] memory[i];
	}
	delete[] memory;
	rows = 0; cols = nullptr;
	memory = nullptr;
}

template <class type>
matrix<type>::matrix(){
	rows = 0; cols = 0;
	memory = nullptr;
}
template <class type>
matrix<type>::matrix(int x, int y){
	int i;
	rows = x;
	cols = y;
	if(rows==0||cols==0){
		memory = nullptr;
	}
	else{
		memory = new type*[rows];
		for(i=0;i<rows;i++){
			memory[i] = new type[cols];
		}
	}
}
template <class type>
matrix<type>::matrix(int x, int y, type z){
	int i, j;
	rows = x;
	cols = y;
	if(rows==0||cols==0){
		memory = nullptr;
	}
	else{
		memory = new type*[rows];
		for(i=0;i<rows;i++){
			memory[i] = new type[cols];
		}
		for(i=0;i<rows;i++){
			for(j=0;j<cols;j++){
				memory[i][j] = z;
			}
		}
	}
}
template <class type>
matrix<type>::~matrix(){
	int i;
	for(i=0;i<rows;i++){
		delete[] memory[i];
	}
	delete[] memory;
	rows = 0; cols = 0;
	memory = nullptr;
}

template <class type>
array<type>::array(){
	length = 0;
	memory = nullptr;
}
template <class type>
array<type>::array(int x){
	length = x;
	if(length==0){
		memory = nullptr;
	}
	else{
		memory = new type[length];
	}
}
template <class type>
array<type>::array(int x, type z){
	int i;
	length = x;
	if(length==0){
		memory = nullptr;
	}
	else{
		memory = new type[length];
		for(i=0;i<length;i++){
			memory[i] = z;
		}
	}
}
template <class type>
array<type>::~array(){
	delete[] memory;
	length = 0;
	memory = nullptr;
}

/*********************
	Retrieve data
*********************/
template <class type>
type matrix4D<type>::value(int x, int y, int z, int d){
	return memory[x][y][z][d];
}
template <class type>
inline void matrix4D<type>::print(int x, int y, int z, int d){
	std::cout << memory[x][y][z][d];
}

template <class type>
type storage<type>::value(int x, int y){
	return memory[x][y];
}
template <class type>
inline void storage<type>::print(int x, int y){
	std::cout << memory[x][y];
}

template <class type>
type matrix<type>::value(int x, int y){
	return memory[x][y];
}
template <class type>
inline void matrix<type>::print(int x, int y){
	std::cout << memory[x][y];
}


template <class type>
type array<type>::value(int x){
	return memory[x];
}
template <class type>
inline void array<type>::print(int x){
	std::cout << memory[x];
}
/*********************
	Alter data
*********************/
template <class type>
inline void matrix4D<type>::set(int x, int y, int z, int d, type a){
	memory[x][y][z][d] = a;
}
template <class type>
void matrix4D<type>::resize(int x, int y, int z, int d, bool a){
	int i, j, k, l;
	if(a==true){
		type**** temp = new type***[x];
		for(i=0;i<x;i++){
			temp[i] = new type**[y];
			for(j=0;j<y;j++){
				temp[i][j] = new type*[z];
				for(k=0;k<z;k++){
					temp[i][j][k] = new type[d];
				}
			}
		}
		for(i=0;i<(x<d1?x:d1);i++){
			for(j=0;j<(y<d2?y:d2);j++){
				for(k=0;k<(z<d3?z:d3);k++){
					for(l=0;l<(d<d4?d:d4);l++){
						temp[i][j] = memory[i][j];
					}
				}
			}
		}

		for(i=0;i<d1;i++){
			for(j=0;j<d2;j++){
				for(k=0;k<d3;k++){
					delete[] memory[i][j][k];
				}
				delete[] memory[i][j];
			}
			delete[] memory[i];
		}
		delete[] memory;
		d1 = x;	d2 = y;
		d3 = z; d4 = d;
		memory = new type***[d1];
		for(i=0;i<d1;i++){
			memory[i] = new type**[d2];
			for(j=0;j<d2;j++){
				memory[i][j] = new type*[d3];
				for(k=0;k<d3;k++){
					memory[i][j][k] = new type[d4];
				}
			}
		}

		for(i=0;i<d1;i++){
			for(j=0;j<d2;j++){
				for(k=0;k<d3;k++){
					for(l=0;l<d4;l++){
						memory[i][j] = temp[i][j];
					}
				}
			}
		}

		for(i=0;i<d1;i++){
			for(j=0;j<d2;j++){
				for(k=0;k<d3;k++){
					delete[] temp[i][j][k];
				}
				delete[] temp[i][j];
			}
			delete[] temp[i];
		}
		delete[] temp;
	}
	else if(a==false){
		for(i=0;i<d1;i++){
			for(j=0;j<d2;j++){
				for(k=0;k<d3;k++){
					delete[] memory[i][j][k];
				}
				delete[] memory[i][j];
			}
			delete[] memory[i];
		}
		delete[] memory;
		d1 = x;	d2 = y;
		d3 = z; d4 = d;
		memory = new type***[d1];
		for(i=0;i<d1;i++){
			memory[i] = new type**[d2];
			for(j=0;j<d2;j++){
				memory[i][j] = new type*[d3];
				for(k=0;k<d3;k++){
					memory[i][j][k] = new type[d4];
				}
			}
		}
	}
}

template <class type>
inline void storage<type>::set(int x, int y, type z){
	memory[x][y] = z;
}
template <class type>
void storage<type>::resize(int x, int*y, bool z){
	int i, j;
	if(z==true){
		type** temp = new type*[x];
		for(i=0;i<x;i++){
			temp[i] = new type[y[i]];
		}
		for(i=0;i<(x<rows?x:rows);i++){
			for(j=0;j<(y[i]<cols[i]?y[i]:cols[i]);j++){
				temp[i][j] = memory[i][j];
			}
		}

		for(i=0;i<rows;i++){
			delete[] memory[i];
		}
		delete[] memory;
		delete[] cols;

		rows = x;
		cols = new int[rows];
		for(i=0;i<rows;i++){
			cols[i] = y[i];
		}
		memory = new type*[rows];
		for(i=0;i<rows;i++){
			memory[i] = new type[cols[i]];
		}

		for(i=0;i<rows;i++){
			for(j=0;j<cols[i];j++){
				memory[i][j] = temp[i][j];
			}
		}

		for(i=0;i<rows;i++){
			delete[] temp[i];
		}
		delete[] temp;
	}
	else if(z==false){
		for(i=0;i<rows;i++){
			delete[] memory[i];
		}
		delete[] memory;
		delete[] cols;

		rows = x;
		cols = new int[rows];
		for(i=0;i<rows;i++){
			cols[i] = y[i];
		}
		memory = new type*[rows];
		for(i=0;i<rows;i++){
			memory[i] = new type[cols[i]];
		}
	}
}

template <class type>
inline void matrix<type>::set(int x, int y, type z){
	memory[x][y] = z;
}
template <class type>
void matrix<type>::resize(int x, int y, bool z){
	int i, j;
	if(z==true){
		type** temp = new type*[x];
		for(i=0;i<x;i++){
			temp[i] = new type[y];
		}
		for(i=0;i<(x<rows?x:rows);i++){
			for(j=0;j<(y<cols?y:cols);j++){
				temp[i][j] = memory[i][j];
			}
		}

		for(i=0;i<rows;i++){
			delete[] memory[i];
		}
		delete[] memory;
		rows = x;
		cols = y;
		memory = new type*[rows];
		for(i=0;i<rows;i++){
			memory[i] = new type[cols];
		}

		for(i=0;i<rows;i++){
			for(j=0;j<cols;j++){
				memory[i][j] = temp[i][j];
			}
		}

		for(i=0;i<rows;i++){
			delete[] temp[i];
		}
		delete[] temp;
	}
	else if(z==false){
		for(i=0;i<rows;i++){
			delete[] memory[i];
		}
		delete[] memory;
		rows = x;
		cols = y;
		memory = new type*[rows];
		for(i=0;i<rows;i++){
			memory[i] = new type[cols];
		}
	}
}

template <class type>
inline void array<type>::set(int x, type z){
	memory[x] = z;
}
template <class type>
void array<type>::resize(int x, bool z){
	int i;
	if(z==true){
		type* temp = new type[x];
		for(i=0;i<(x<length?x:length);i++){
			temp[i] = memory[i];
		}

		delete[] memory;
		length = x;
		memory = new type[length];

		for(i=0;i<length;i++){
			memory[i] = temp[i];
		}

		delete[] temp;
	}
	else if(z==false){
		delete[] memory;
		length = x;
		memory = new type[length];
	}
}

#endif
