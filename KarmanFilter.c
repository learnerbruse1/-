#include <stdlib.h>

double ** KarmanFilter(int row, int col, double (*s)[col]);       //在导入的数组s中 row == 行,col == 列
double ** zeros(int row, int col);                              //创建单位矩阵
double ** Transpo(int row, int col, double (*s)[col]);          //转置函数
void Free_Mat(double ** t, int row);                               //释放内存
int length(int row, int col);
double ** Mat_Mul(int row1, int col1, double ** arr1, int col2, double ** arr2); //矩阵相乘第一个的列数等于第二个的行数(col1 == row2)
double ** Mat_Mul_11(int row1, int col1, double (*arr1)[col1], int col2, double (*arr2)[col2]);
double ** Mat_Mul_12(int row1, int col1, double (*arr1)[col1], int col2, double ** arr2);
double ** Mat_Mul_21(int row1, int col1, double ** arr1, int col2, double (*arr2)[col2]);
double ** EYE(int m, int n);      //生成单位矩阵

double ** KarmanFilter(int row, int col, double (*s)[col]){      // col == 1
    double ** st;                        //转置后的数组
    double ** Y;                            
    double ** B;
    double ** At;
    double ** Ht;
    double ** C;

    st = Transpo(row, col, s);
    double T = (double)0.1;
    int m = length(row, col);
    Y = zeros(2, m+1);
    double Y0[2][1] = {{0},
                      {1}
    };
    for (int i = 0; i < 2; i++){
        Y[i][0] = Y0[i][0];
    }
    double A[2][2] = {{1, T},
                     {0, 1}
    };
    At = Transpo(2, 2, A);
    double Bt[1][2] = {0.5 * (T * T), T
    };
    B = Transpo(1, 2, Bt);
    double H[1][2] = {1, 0
    };
    Ht = Transpo(1, 2, H);
    double C0[2][2] = {{0, 0},
                      {0, 1}
    };
    C = zeros(2, 2*m + 2);                   //= {C0, zeros(2,2*m)
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
                C[i][j] = C0[i][j];
        }
    }
    double Q = (0.25) * (0.25);
    double R = (0.25) * (0.25);

    double C_Change[2][2];
    double ** K_Gain1;
    double ** K_Gain2;
    double ** K_Gain22;
    double ** K;
    double Y_Change[2][1];
    double ** Y_Cal2;
    double ** Y_Cal22;
    double ** Y_For;
    double ** eye;
    double ** eyee;
    double ** C_Var;
    double ** C_Con1;
    double ** C_Con11;
    double ** C_Con2;
    for (int a = 0; a < 2; a++)
        B[a][0] *= Q; 
    C_Con2 = Mat_Mul_21(2, 1, B, 2, Bt);
    
    for (int n = 0; n < m; n++){        //C语言的数组下标从0开始
        
        int i = n * 2;

        //卡尔曼增益矩阵
        for (int a = 0; a < 2; a++){
                C_Change[a][0] = C[a][i];
                C_Change[a][1] = C[a][i+1];
        }
        K_Gain1 = Mat_Mul_12(2, 2, C_Change, 1, Ht);
        K_Gain2 = Mat_Mul_11(1, 2, H, 2, C_Change);
        K_Gain22 = Mat_Mul(1, 2, K_Gain2, 1, Ht);
        K_Gain22[0][0] += R;
        K_Gain22[0][0] = 1.0 / K_Gain22[0][0];    //该矩阵仅为一阶，所以逆矩阵为该元素的倒数
        K = Mat_Mul(2, 1, K_Gain1, 1, K_Gain22); 

        //状态估计                   
        for (int a = 0; a < 2; a++)
            Y_Change[a][0] = Y[a][n];
        Y_Cal2 = Mat_Mul_11(1, 2, H, 1, Y_Change);
        Y_Cal2[0][0] = st[0][n] - Y_Cal2[0][0];     //Y_Cal2和st仅有一行
        Y_Cal22 = Mat_Mul(2, 1, K, 1, Y_Cal2);
        for(int a = 0; a < 2; a++){
            Y[a][n] += Y_Cal22[a][0];
        }
        for (int a = 0; a < 2; a++)
            Y_Change[a][0] = Y[a][n];
        
        //状态预测
        Y_For = Mat_Mul_11(2, 2, A, 1, Y_Change);
        for (int a = 0; a < 2; a++){ 
            Y[a][n+1] = Y_For[a][0];

        }
        
        //估计方差阵
        eyee = Mat_Mul_21(2, 1, K, 2, H);
        eye = EYE(2, 2);
        for (int i = 0; i < 2; i++){
            for (int j = 0; j < 2; j++)
                eyee[i][j] = eye[i][j] - eyee[i][j];
        }
        C_Var = Mat_Mul_21(2, 2, eyee, 2, C_Change);
        for (int a = 0; a < 2; a++){
            C[a][i] = C_Var[a][0];
            C[a][i+1] = C_Var[a][1];
        }
        for (int a = 0; a < 2; a++){
                C_Change[a][0] = C[a][i];
                C_Change[a][1] = C[a][i+1];
        }

        //进一步预测方差阵
        C_Con1 = Mat_Mul_11(2, 2, A, 2, C_Change);
        C_Con11 = Mat_Mul(2, 2, C_Con1, 2, At);
        
        for (int a = 0; a < 2; a++){
            C[a][i+2] = C_Con11[a][0] + C_Con2[a][0];
            C[a][i+3] = C_Con11[a][1] + C_Con2[a][1];
            
        }
        
        Free_Mat(K_Gain1, 2);
        Free_Mat(K_Gain2, 1);
        Free_Mat(K_Gain22, 1);
        Free_Mat(K,2);
        Free_Mat(Y_Cal2, 1);
        Free_Mat(Y_Cal22, 2);
        Free_Mat(Y_For, 2);
        Free_Mat(eye, 2);
        Free_Mat(eyee, 2);
        Free_Mat(C_Var, 2);
        Free_Mat(C_Con1, 2);
        Free_Mat(C_Con11, 2);
        
    }
    
    double ** yp = zeros(1, m+1);
    for (int i = 0; i <= m; i++){
        yp[0][i] = Y[0][i];
    }

    Free_Mat(st, col);
    Free_Mat(C_Con2, 2);
    Free_Mat(Y,2);
    Free_Mat(B, 2);
    Free_Mat(At, 2);
    Free_Mat(Ht, 2);
    Free_Mat(C, 2);

    return yp;            //分配了动态内存
}

double ** Transpo(int row, int col, double (*s)[col]){              //对数组进行转置操作
    double ** t = (double**) malloc(col * sizeof(double*));
    for(int i = 0; i < col; i++){
        t[i] = (double*) malloc(row * sizeof(double));
    }
    for(int i = 0; i < row; i++){
        for (int j = 0; j < col; j++){
            t[j][i] = s[i][j];
        }
    }
    return t;
}

void Free_Mat(double ** t, int row){                                 //释放数组的动态内存
    for (int i = 0; i < row; i++)
        free(t[i]);
    free(t);
}

int length(int row, int col){                                       //返回数组较长的 列/行 数
    if (row < col)
        return col;
    else
        return row;
}

double ** zeros(int row, int col){                                  //创建存储数据全为0的数组
    double ** arr = (double**) calloc(row , sizeof(double*)); //进行内存分配
    for(int i = 0; i < row; i++){
        arr[i] = (double*) calloc(col , sizeof(double));
    }
    return arr;
}

double ** Mat_Mul(int row1, int col1, double ** arr1, int col2, double ** arr2){    //二维指针的数组与二维指针的数组相乘
    double ** result = zeros(row1,col2);                        //调用zeros函数      //矩阵相乘第一个矩阵的列数必须等于第二个矩阵的行数  即col1 = row2
    for(int i = 0; i < row1; i++){                                                  //结果的 行数=row1    列数=col2  
        for(int j = 0; j < col2; j++){                                                  
            for (int k = 0; k < col1; k++){
                result[i][j] += arr1[i][k] * arr2[k][j];        //计算矩阵相乘结果
            }
        }
    }
    return result;
}

double ** Mat_Mul_11(int row1, int col1, double (*arr1)[col1], int col2, double (*arr2)[col2]){  //数组指针(存储方式的指针只有一维[连续存储])的数组与数组指针的数组相乘 
    double ** result = zeros(row1,col2);        
    for(int i = 0; i < row1; i++ ){
        for(int j = 0; j < col2; j++){
            for (int k = 0; k < col1; k++){
                result[i][j]  += arr1[i][k] * arr2[k][j];  
            }
        }
    }
    return result;
}

double ** Mat_Mul_12(int row1, int col1, double (*arr1)[col1], int col2, double ** arr2){   //数组指针的数组与二维指针的数组相乘
    double ** result = zeros(row1,col2);
    for(int i = 0; i < row1; i++ ){
        for(int j = 0; j < col2; j++){
            for (int k = 0; k < col1; k++){
                result[i][j]  += arr1[i][k] * arr2[k][j];
            }
        }
    }
    return result;
}

double ** Mat_Mul_21(int row1, int col1, double ** arr1, int col2, double (*arr2)[col2]){   //二维指针的数组与数组指针的数组相乘
    double ** result = zeros(row1,col2);
    for(int i = 0; i < row1; i++ ){
        for(int j = 0; j < col2; j++){
            for (int k = 0; k < col1; k++){
                result[i][j]  += arr1[i][k] * arr2[k][j];
            }
        }
    }
    return result;
}

double ** EYE(int m, int n){            //创建单位矩阵
    double ** Mat = zeros(m, n);
    int i;
    if (m >= n)                 //比较行数和列数大小
        i = n;
    else
        i = m;
    for (int a = 0; a < i; a++)
        Mat[a][a] = (double)1;
    return Mat;
}

//需对返回的指针free
/*
void Free_Mat(float **yp, int row){         // row == 1
    for (int i = 0; i < row; i++)
        free(yp[i]);
    free(yp);
}
*/