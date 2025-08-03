// MAKE THIS INLINE.

#ifndef LINALG_HPP
#define LINALG_HPP

#include "Constants/Constants.hpp"
#include <cmath>
#include <Accelerate/Accelerate.h>

namespace Planck {
    inline complex* getRawArray(vector& a){
        return a.data();
    }

    inline vector fromRaw(complex* aRaw, int n){
        vector a; 
        for(int i = 0; i < n; i++) a.push_back(aRaw[i]);

        return a;
    }

    inline matrix fromRaw(complex* ARaw, int n, int m){
        matrix A(n);
        for(int i = 0; i < n; i++) 
            for(int j = 0; j < m; j++)
                A[i].push_back(ARaw[i * m + j]);

        return A;
    }

    inline matrix zeros(int n, int m) {
        matrix A(n);
        for(int i = 0; i < n; i++) A[i].resize(m);

        return A; 
    }

    inline vector zeros(int n) {
        return vector(n);
    }

    inline vector ones(int n){
        vector a = zeros(n);
        for(int i = 0; i < n; i++) a[i] = 1; 
        
        return a;
    }

    inline matrix eye(int n){
        matrix A = zeros(n,n);
        for(int i = 0; i < n; i++){
            A[i][i] = 1; 
        }
        return A;
    }

    inline vector flatten(matrix A){
        vector flattened;
        for(int i = 0; i < A.size(); i++){
            flattened.insert(flattened.end(), A[i].begin(), A[i].end());
        }

        return flattened;
    }

    inline matrix reshape(vector a, int n, int m){
        matrix A = zeros(n, m);

        for(int i = 0; i < n; i++)
            for(int j = 0; j < m; j++)
                A[i][j] = a[i * m + j];

        return A;
    }

    inline vector operator+(vector a, vector b){
        complex* aRaw = getRawArray(a);
        complex* bRaw = getRawArray(b);

        complex alpha = {1.0, 0.0};

        cblas_zaxpy(a.size(), &alpha, aRaw, 
                1, bRaw, 1);

        return fromRaw(bRaw, a.size());
    }

    inline vector operator-(vector a, vector b){
        complex* aRaw = getRawArray(a);
        complex* bRaw = getRawArray(b);

        complex alpha = {-1.0, 0.0};

        cblas_zaxpy(a.size(), &alpha, bRaw, 
                1, aRaw, 1);

        return fromRaw(aRaw, a.size());
    }

    inline vector operator*(vector a, complex b){
        complex* aRaw = getRawArray(a);

        cblas_zscal(a.size(), &b, aRaw, 1);

        return fromRaw(aRaw, a.size());
    }

    inline vector operator*(complex a, vector b){
        return b * a;
    }

    inline vector operator/(vector a, complex b){
        complex* aRaw = getRawArray(a);

        complex c = complex(1)/b;
        cblas_zscal(a.size(), &c, aRaw, 1);

        return fromRaw(aRaw, a.size());
    }

    inline vector operator/(complex a, vector b){
        return b / a;
    }

    inline matrix operator+(matrix A, matrix B){
        vector AFlat = flatten(A);
        vector BFlat = flatten(B);

        return reshape(AFlat + BFlat, A.size(), A[0].size());
    }

    inline matrix operator-(matrix A, matrix B){
        vector AFlat = flatten(A);
        vector BFlat = flatten(B);

        return reshape(AFlat - BFlat, A.size(), A[0].size());
    }

    inline matrix operator*(matrix A, complex b){
        vector AFlat = flatten(A);

        return reshape(AFlat * b, A.size(), A[0].size());
    }

    inline matrix operator*(complex a, matrix B){
        vector BFlat = flatten(B);

        return reshape(a * BFlat, B.size(), B[0].size());
    }

    inline matrix operator/(matrix A, complex b){
        vector AFlat = flatten(A);

        return reshape(AFlat / b, A.size(), A[0].size());
    }

    inline matrix operator/(complex a, matrix B){
        vector BFlat = flatten(B);

        return reshape(BFlat / a, B.size(), B[0].size());
    }

    inline vector operator+(vector a){
        return a; 
    }

    inline vector operator-(vector b){
        return -1 * b; 
    }

    inline matrix operator+(matrix A){
        return A; 
    }

    inline matrix operator-(matrix B){
        return -1 * B;
    }

    inline void operator+=(matrix& A, matrix B){
        A = A + B; 
    }

    inline void operator-=(matrix& A, matrix B){
        A = A - B; 
    }

    inline matrix dot(matrix A, matrix B){
        vector AFlat = flatten(A);
        vector BFlat = flatten(B);

        complex* ARaw = getRawArray(AFlat);
        complex* BRaw = getRawArray(BFlat);

        vector C = zeros(A.size() * B[0].size());
        complex* CRaw = getRawArray(C);

        complex alpha = {1.0, 0.0};
        complex beta = {0.0, 0.0};

        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
            A.size(), B[0].size(), B.size(), &alpha,
            ARaw, A[0].size(), BRaw, B[0].size(), 
            &beta, CRaw, B[0].size());

        return fromRaw(CRaw, A.size(), B[0].size());
    }

    inline vector dot(matrix A, vector b){

        vector AFlat = flatten(A);

        complex* ARaw = getRawArray(AFlat);
        complex* bRaw = getRawArray(b);

        vector c = zeros(A.size());
        complex* cRaw = getRawArray(c);

        complex alpha = {1.0, 0.0};
        complex beta = {0.0, 0.0};

        cblas_zgemv(CblasRowMajor, CblasNoTrans, 
            A.size(), A[0].size(), &alpha,
            ARaw, A[0].size(), bRaw, 1, 
            &beta, cRaw, 1);

        return fromRaw(cRaw, A.size());

    }

    inline vector dot(vector a, matrix B){
        vector BFlat = flatten(B);

        complex* aRaw = getRawArray(a);
        complex* BRaw = getRawArray(BFlat);

        vector c = zeros(B[0].size());
        complex* cRaw = getRawArray(c);

        complex alpha = {1.0, 0.0};
        complex beta = {0.0, 0.0};

        cblas_zgemv(CblasRowMajor, CblasTrans, 
            B.size(), B[0].size(), &alpha,
            BRaw, B[0].size(), aRaw, 1, 
            &beta, cRaw, 1);

        return fromRaw(cRaw, B[0].size());
    }

    inline complex dot(vector a, vector b){
        complex* aRaw = getRawArray(a);
        complex* bRaw = getRawArray(b);

        complex cRaw;

        cblas_zdotu_sub(a.size(), aRaw, 1, 
            bRaw, 1, &cRaw);

        return cRaw;
    }

    inline vector kron(vector a, vector b){
        complex* aRaw = getRawArray(a);
        complex* bRaw = getRawArray(b);

        vector c = zeros(a.size() * b.size());
        complex* cRaw = getRawArray(c);

        for(int i = 0; i < a.size(); i++){
            vector temp = b; 
            complex* tempRaw = getRawArray(temp);

            cblas_zscal(temp.size(), &aRaw[i], tempRaw, 1);

            cblas_zcopy(temp.size(), tempRaw, 1, 
                cRaw + temp.size() * i, 1);
        }
        return fromRaw(cRaw, c.size());
    }

    inline matrix kron(matrix A, matrix B){
        vector AFlat = flatten(A);

        vector BFlat = flatten(B);

        complex* ARaw = getRawArray(AFlat);
        complex* BRaw = getRawArray(BFlat);

        vector C = zeros(A.size() * B.size() * A[0].size() * B[0].size());
        complex* CRaw = getRawArray(C);

        for(int i = 0; i < A.size(); i++){
            for(int j = 0; j < A[0].size(); j++){
                vector temp = BFlat; 
                complex* tempRaw = getRawArray(temp);

                cblas_zscal(temp.size(), &ARaw[i * A[0].size() + j], tempRaw, 1);

                // move this into the correct position.
                for(int k = 0; k < B.size(); k++){
                    cblas_zcopy(B[0].size(), tempRaw + k * B[0].size(), 1, 
                        CRaw + i * A[0].size() * B.size() * B[0].size() + j * B[0].size() + k * A[0].size() * B[0].size(), 1);
                }
            }
        }

        return fromRaw(CRaw, A.size() * B.size(), A[0].size() * B[0].size());
    }

    inline vector hermitian(vector a){
        for(int i = 0; i < a.size(); i++) a[i] = conj(a[i]);
        return a; 
    }

    inline vector linspace(double start, double end, int num){
        double delta = (end - start) / (double(num) - 1);

        vector a;
        for(int i = start; i <= end; i+=delta){
            a.push_back(i);
        }

        return a;
    }

    inline double norm(vector a) {
        return std::sqrt(dot(a, a)).real();
    }

    inline double norm(matrix A){
        vector AFlat = flatten(A);

        return norm(AFlat);
    }
}





#endif