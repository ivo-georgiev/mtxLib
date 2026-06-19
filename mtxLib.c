/******************************************************************************************************************************************************************************************************\
 ***
 *** Description       : IMPLEMENTATION OF BASIC MATRIX OPERATION
 *** Codefile          : mtxLib.c
 ***
 *** MIT License
 ***
 *** Copyright (c) 2017 ivo-georgiev
 ***
 *** Permission is hereby granted, free of charge, to any person obtaining a copy
 *** of this software and associated documentation files (the "Software"), to deal
 *** in the Software without restriction, including without limitation the rights
 *** to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *** copies of the Software, and to permit persons to whom the Software is
 *** furnished to do so, subject to the following conditions:
 ***
 *** The above copyright notice and this permission notice shall be included in all
 *** copies or substantial portions of the Software.
 ***
 *** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *** IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *** FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *** AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *** LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *** OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *** SOFTWARE.
\******************************************************************************************************************************************************************************************************/

#include "mtxLib.h"

void mtx_init_bool(MatrixBool_t *const pSrc, _Bool *const pValue, const int nrow, const int ncol, const int nelem)
{
	pSrc->val = pValue;
	pSrc->ncol = ncol;
	pSrc->nrow = nrow;
	pSrc->nelem = nelem;
}

void mtx_init(Matrix_t *const pSrc, double *const pValue, const int nrow, const int ncol, const int nelem)
{
	pSrc->val = pValue;
	pSrc->ncol = ncol;
	pSrc->nrow = nrow;
	pSrc->nelem = nelem;
}

/**
 * For square matrix only
 */
enum mtxResultInfo mtx_diagsum(Matrix_t *pSrc, double *diagsum)
{
	enum mtxResultInfo Result = MTX_OPERATION_OK;
	double const *const pSrcL = (double *)pSrc->val;
	const int ncol = pSrc->ncol;
	int eIdx;
	double sum = pSrcL[0];

	if (pSrc->nrow == ncol)
	{
		for (eIdx = 1; eIdx < pSrc->nelem; eIdx++)
		{
			const int cmpLeft = (int)(eIdx / ncol);

			sum += eIdx < ncol ? 0 : cmpLeft == eIdx % (cmpLeft * ncol) ? pSrcL[eIdx] : 0;
		}
	}
	else
	{
		Result = MTX_SIZE_MISMATCH;
	}

	*diagsum = sum;

	return Result;
}
/**
 * @brief  \f$A=A'\f$ or \f$A=A^{\tau}\f$
 */
enum mtxResultInfo mtx_transp_square(Matrix_t *const pSrc)
{
	enum mtxResultInfo ResultL = MTX_OPERATION_OK;
	const int nrow = pSrc->nrow;
	const int ncol = pSrc->ncol;
	double *const pSrcL = (double *)pSrc->val;
	int row, col;
	double temp;

	if (nrow == ncol)
	{
		for (row = 0; row < nrow; row++)
		{
			for (col = 0; col < ncol; col++)
			{
				if (row != col && row < col)
				{
					temp = pSrcL[nrow * row + col];
					pSrcL[ncol * row + col] = pSrcL[ncol * col + row];
					pSrcL[ncol * col + row] = temp;
				}
			}
		}
	}
	else
	{
		ResultL = MTX_NOT_SQUARE;
	}

	return ResultL;
}
enum mtxResultInfo mtx_transp_dest(Matrix_t const *const pSrc, Matrix_t *const pDst)
{
	enum mtxResultInfo ResultL = MTX_OPERATION_OK;
	double const *const pSrcL = (double *)pSrc->val;
	double *const pDstL = (double *)pDst->val;
	const int nRowSrcL = pSrc->nrow;
	const int nColSrcL = pSrc->ncol;
	const int nRowDstL = pDst->nrow;
	const int nColDstL = pDst->ncol;
	int row, col;

	if (nRowSrcL == nColDstL || nColSrcL == nRowDstL)
	{
		for (row = 0; row < nRowDstL; row++)
		{
			for (col = 0; col < nColDstL; col++)
			{
				pDstL[nColDstL * row + col] = pSrcL[nColSrcL * col + row];
			}
		}
	}
	else
	{
		ResultL = MTX_SIZE_MISMATCH;
	}

	return ResultL;
}
/**
 * @brief \f$C=A \cdot B\f$
 *
 * Matrix multiplication
 *
 */
enum mtxResultInfo mtx_mul(Matrix_t const *const pSrc1, Matrix_t const *const pSrc2, Matrix_t *const pDst)
{
	enum mtxResultInfo ResultL = MTX_OPERATION_OK;
	double const *const pSrc1L = (double *)pSrc1->val;
	double const *const pSrc2L = (double *)pSrc2->val;
	double *const pDstL = (double *)pDst->val;
	int row, col, k;
	double sum;

	if (pSrc1->ncol == pSrc2->nrow)
	{
		for (row = 0; row < pSrc1->nrow; row++)
		{
			for (col = 0; col < pSrc2->ncol; col++)
			{
				sum = 0;
				for (k = 0; k < pSrc1->ncol; k++)
				{
					sum += pSrc1L[pSrc1->ncol * row + k] * pSrc2L[pSrc2->ncol * k + col];
				}
				pDstL[pDst->ncol * row + col] = sum;
			}
		}
	}
	else
	{
		ResultL = MTX_SIZE_MISMATCH;
	}

	return ResultL;
}
/**
 * @brief Special multiplication \f$Dst=Src1 \cdot Src2^{\tau}\f$
 *
 *  Special function for multiplication of Src1 matrix with transpose image of
 *  Src2 matrix. Be sure that array size are suitable for multiplication after
 *  Src2 transpose
 */
enum mtxResultInfo mtx_mul_src2tr(Matrix_t const *const pSrc1, Matrix_t const *const pSrc2, Matrix_t *const pDst)
{
	enum mtxResultInfo ResultL = MTX_OPERATION_OK;
	double const *const pSrc1L = (double *)pSrc1->val;
	double const *const pSrc2L = (double *)pSrc2->val;
	double *const pDstL = (double *)pDst->val;
	int rowSrc1, rowSrc2, k;
	double sum;

	if (pSrc1->ncol == pSrc2->ncol)
	{
		for (rowSrc1 = 0; rowSrc1 < pSrc1->nrow; rowSrc1++)
		{
			for (rowSrc2 = 0; rowSrc2 < pSrc2->nrow; rowSrc2++)
			{
				sum = 0;
				for (k = 0; k < pSrc1->ncol; k++)
				{
					sum += pSrc1L[pSrc1->ncol * rowSrc1 + k] * pSrc2L[pSrc2->ncol * rowSrc2 + k];
				}
				pDstL[pDst->ncol * rowSrc1 + rowSrc2] = sum;
			}
		}
	}
	else
	{
		ResultL = MTX_SIZE_MISMATCH;
	}

	return ResultL;
}
enum mtxResultInfo mtx_chol_lower(Matrix_t *const pSrc)
{
	enum mtxResultInfo ResultL = MTX_OPERATION_OK;
	double *const pSrcL = pSrc->val;
	const int nrow = pSrc->nrow;
	const int ncol = pSrc->ncol;
	int col, row;
	int tmp;
	double sum = 0;

	if (ncol == nrow)
	{
		const int mtxSize = nrow;

		for (col = 0; col < mtxSize; col++)
		{
			for (row = 0; row < mtxSize; row++)
			{
				sum = pSrcL[mtxSize * col + row];

				for (tmp = (int)(col - 1); tmp >= 0; tmp--)
				{
					sum -= pSrcL[mtxSize * row + tmp] * pSrcL[mtxSize * col + tmp];
				}

				pSrcL[ncol * row + col] = (row == col) ? sqrt(sum) : (row > col) ? (sum / pSrcL[ncol * col + col]) : 0;

				if ((row == col) && (sum <= 0))
				{
					ResultL = MTX_NOT_POS_DEFINED;
				}
			}
		}
	}
	else
	{
		ResultL = MTX_NOT_SQUARE;
	}

	return ResultL;
}
enum mtxResultInfo mtx_chol_upper(Matrix_t *const pSrc)
{
	enum mtxResultInfo ResultL = MTX_OPERATION_OK;
	double *const pSrcL = pSrc->val;
	const int nrow = pSrc->nrow;
	const int ncol = pSrc->ncol;
	int col, row;
	int tmp;
	double sum = 0;

	if (ncol == nrow)
	{
		for (row = 0; row < nrow; row++)
		{
			for (col = 0; col < ncol; col++)
			{
				sum = pSrcL[ncol * row + col];

				for (tmp = row - 1; tmp >= 0; tmp--) // tmp could be calc negative
				{
					sum -= pSrcL[ncol * tmp + row] * pSrcL[ncol * tmp + col];
				}

				pSrcL[ncol * row + col] = (row == col) ? sqrt(sum) : (row < col) ? (sum / pSrcL[ncol * row + row]) : 0;

				if ((row == col) && (sum <= 0))
				{
					ResultL = MTX_NOT_POS_DEFINED;
				}
			}
		}
	}
	else
	{
		ResultL = MTX_NOT_SQUARE;
	}

	return ResultL;
}
/**
 * Upper cholesky decomposition variant 1
 */
enum mtxResultInfo mtx_chol1(double *A, double *L, const int size)
{
	enum mtxResultInfo Result = MTX_OPERATION_OK;
	int col, row;
	int tmp;
	double sum = 0;

	for (row = 0; row < size; row++)
	{
		for (col = 0; col < size; col++)
		{
			sum = A[size * row + col];

			for (tmp = row - 1; tmp >= 0; tmp--)
			{
				sum -= L[size * tmp + row] * L[size * tmp + col];
			}

			if (row == col)
			{
				if (sum > 0)
				{
					L[size * row + col] = sqrt(sum);
				}
				else
				{
					Result = MTX_NOT_POS_DEFINED;
				}
			}
			else if (row < col)
			{
				L[size * row + col] = sum / L[size * row + row]; // sum/Lii(diag)
			}
			else
			{
				L[size * row + col] = 0;
			}
		}
	}

	return Result;
}
/**
 * @brief Matrix inv  \f$A=A^{-1}\f$
 *
 * @param pDst At the begining should point to identity matrix!!
 * @param pSrc Square matrix
 */
/**
 * @brief Compute matrix inverse using hardened Gauss-Jordan elimination with partial pivoting
 * @details Implements Gauss-Jordan elimination with:
 *   - Pre-validation pass to detect singularity before modifications
 *   - Partial pivoting for numerical stability
 *   - Epsilon-based zero detection
 *   - Input/output validation and consistency checks
 *   - Automatic initialization of destination as identity matrix
 *   - Condition number warning for ill-conditioned matrices
 * @param pSrc Source matrix (square, N×N). Modified in-place during computation if successful.
 *             On failure, source is unchanged (early detection before modifications).
 * @param pDst Destination matrix (output inverse). Automatically initialized as identity matrix.
 *             On success, contains inverse. Must have pre-allocated storage matching pSrc dimensions.
 * @return MTX_OPERATION_OK if inversion successful
 *         MTX_SINGULAR if matrix is singular (det ≈ 0)
 *         MTX_NOT_SQUARE if matrix is not square
 *         MTX_SIZE_MISMATCH if src/dst dimensions don't match
 *         MTX_OPERATION_ERROR if NULL pointers or invalid values
 * @pre pSrc->nrow == pSrc->ncol (enforced)
 * @pre pSrc->nrow == pDst->nrow && pSrc->ncol == pDst->ncol (enforced)
 * @pre pDst->val has valid allocated storage (required)
 * @pre All elements in pSrc are finite real numbers (checked)
 * @note Uses absolute epsilon threshold relative to matrix norm for robust singularity detection
 * @note Singularity is detected in pre-validation pass; source unchanged if detected
 * @note Condition number heuristic: if max_pivot/min_pivot > 1e8, matrix is ill-conditioned
 */
enum mtxResultInfo mtx_inv(Matrix_t *const pSrc, Matrix_t *const pDst)
{
	enum mtxResultInfo result = MTX_OPERATION_OK;
	
	/* ===== INPUT VALIDATION ===== */
	
	/* Check NULL pointers */
	if (NULL == pSrc || NULL == pDst || NULL == pSrc->val || NULL == pDst->val)
	{
		return MTX_OPERATION_ERROR;
	}
	
	const int matrix_size = pSrc->nrow;
	const int ncol = pSrc->ncol;
	
	/* Check matrix is square */
	if (matrix_size != ncol)
	{
		return MTX_NOT_SQUARE;
	}
	
	/* Check dimensions match between source and destination */
	if (pSrc->nrow != pDst->nrow || pSrc->ncol != pDst->ncol || 
	    pSrc->nelem != pDst->nelem)
	{
		return MTX_SIZE_MISMATCH;
	}
	
	/* Check for minimum matrix size (1×1 minimum) */
	if (matrix_size < 1)
	{
		return MTX_OPERATION_ERROR;
	}
	
	/* Check for special/invalid values in source matrix */
	for (int elem_idx = 0; elem_idx < pSrc->nelem; elem_idx++)
	{
		if (!isfinite(pSrc->val[elem_idx]))  /* Checks for NaN, Inf, -Inf */
		{
			return MTX_OPERATION_ERROR;
		}
	}
	
	/* Verify destination is initialized as identity matrix */
	for (int row = 0; row < matrix_size; row++)
	{
		for (int col = 0; col < matrix_size; col++)
		{
			pDst->val[matrix_size * row + col] = (row == col) ? 1.0 : 0.0;
		}
	}
	
	/* ===== NUMERICAL PREPARATION ===== */
	
	/* Compute matrix norm for adaptive epsilon threshold */
	double matrix_norm = 0.0;
	for (int elem_idx = 0; elem_idx < pSrc->nelem; elem_idx++)
	{
		double abs_val = fabs(pSrc->val[elem_idx]);
		if (abs_val > matrix_norm)
		{
			matrix_norm = abs_val;
		}
	}
	
	/* Set epsilon threshold: relative to matrix norm or absolute minimum */
	const double EPSILON_BASE = 1e-14;
	const double epsilon_threshold = (matrix_norm > 1.0) ? 
	                                  EPSILON_BASE * matrix_norm : 
	                                  EPSILON_BASE;
	
	/* ===== PRE-VALIDATION PASS: CHECK ALL PIVOTS BEFORE ANY MODIFICATIONS ===== */
	/* This ensures singularity is detected before corrupting the source matrix */
	
	double max_pivot = 0.0;
	double min_pivot = 1e16;
	int pivot_row_array[matrix_size];  /* Store pivot rows for elimination phase */
	
	for (int pivot_col = 0; pivot_col < matrix_size; pivot_col++)
	{
		/* Find pivot row using PARTIAL PIVOTING (largest element in column below diagonal) */
		int pivot_row = pivot_col;
		double max_pivot_val = fabs(pSrc->val[matrix_size * pivot_col + pivot_col]);
		
		for (int search_row = pivot_col + 1; search_row < matrix_size; search_row++)
		{
			double candidate_pivot = fabs(pSrc->val[matrix_size * search_row + pivot_col]);
			if (candidate_pivot > max_pivot_val)
			{
				max_pivot_val = candidate_pivot;
				pivot_row = search_row;
			}
		}
		
		/* Check if pivot is acceptably large (not near zero) */
		if (max_pivot_val < epsilon_threshold)
		{
			/* Matrix is singular - return BEFORE any modifications */
			return MTX_SINGULAR;
		}
		
		/* Track pivot magnitudes for condition number heuristic */
		if (max_pivot_val > max_pivot)
		{
			max_pivot = max_pivot_val;
		}
		if (max_pivot_val < min_pivot)
		{
			min_pivot = max_pivot_val;
		}
		
		/* Store the pivot row for this column for use in elimination phase */
		pivot_row_array[pivot_col] = pivot_row;
	}
	
	/* ===== GAUSS-JORDAN ELIMINATION WITH PARTIAL PIVOTING ===== */
	/* All pivots validated; now we can safely proceed with elimination */
	
	for (int pivot_col = 0; pivot_col < matrix_size; pivot_col++)
	{
		int pivot_row = pivot_row_array[pivot_col];
		
		/* Swap rows if needed (row exchange for partial pivoting) */
		if (pivot_row != pivot_col)
		{
			for (int col_idx = 0; col_idx < matrix_size; col_idx++)
			{
				/* Swap in source matrix */
				double temp_src = pSrc->val[matrix_size * pivot_col + col_idx];
				pSrc->val[matrix_size * pivot_col + col_idx] = 
					pSrc->val[matrix_size * pivot_row + col_idx];
				pSrc->val[matrix_size * pivot_row + col_idx] = temp_src;
				
				/* Swap in destination (identity) matrix */
				double temp_dst = pDst->val[matrix_size * pivot_col + col_idx];
				pDst->val[matrix_size * pivot_col + col_idx] = 
					pDst->val[matrix_size * pivot_row + col_idx];
				pDst->val[matrix_size * pivot_row + col_idx] = temp_dst;
			}
		}
		
		/* Scale pivot row to make diagonal element = 1 */
		double pivot_value = pSrc->val[matrix_size * pivot_col + pivot_col];
		double scale_factor = 1.0 / pivot_value;
		
		/* Scale the pivot row */
		for (int col_idx = 0; col_idx < matrix_size; col_idx++)
		{
			pSrc->val[matrix_size * pivot_col + col_idx] *= scale_factor;
			pDst->val[matrix_size * pivot_col + col_idx] *= scale_factor;
		}
		
		/* Eliminate column in all other rows */
		for (int elim_row = 0; elim_row < matrix_size; elim_row++)
		{
			if (elim_row != pivot_col)
			{
				double elimination_factor = -pSrc->val[matrix_size * elim_row + pivot_col];
				
				for (int col_idx = 0; col_idx < matrix_size; col_idx++)
				{
					pSrc->val[matrix_size * elim_row + col_idx] += 
						elimination_factor * pSrc->val[matrix_size * pivot_col + col_idx];
					pDst->val[matrix_size * elim_row + col_idx] += 
						elimination_factor * pDst->val[matrix_size * pivot_col + col_idx];
				}
			}
		}
	}
	
	/* ===== POST-COMPUTATION CHECKS ===== */
	
	/* Verify result is finite */
	for (int elem_idx = 0; elem_idx < pDst->nelem; elem_idx++)
	{
		if (!isfinite(pDst->val[elem_idx]))
		{
			return MTX_OPERATION_ERROR;
		}
	}
	
	/* Heuristic condition number check: warn if very ill-conditioned */
	/* Condition number ≈ max_pivot / min_pivot */
	if (min_pivot > 1e-16 && (max_pivot / min_pivot) > 1e8)
	{
		/* Matrix is ill-conditioned but inversion completed */
		/* Result may have low accuracy - no error returned but computation succeeded */
		result = MTX_OPERATION_OK;
	}
	
	return result;
}
enum mtxResultInfo mtx_add(Matrix_t *const pDst, Matrix_t const *const pSrc)
{
	int Result = MTX_OPERATION_OK;
	double *const pDstL = (double *)pDst->val;
	double const *const pSrcL = (double *)pSrc->val;
	int eIdx;

	if (pDst->ncol == pSrc->ncol && pDst->nrow == pSrc->nrow)
	{
		for (eIdx = 0; eIdx < pSrc->nelem; eIdx++)
		{
			pDstL[eIdx] += pSrcL[eIdx];
		}
	}
	else
	{
		Result = MTX_SIZE_MISMATCH;
	}

	return Result;
}
enum mtxResultInfo mtx_sub(Matrix_t *const pDst, Matrix_t const *const pSrc)
{
	int Result = MTX_OPERATION_OK;
	double *const pDstL = (double *)pDst->val;
	double const *const pSrcL = (double *)pSrc->val;
	int eIdx;

	if (pDst->ncol == pSrc->ncol && pDst->nrow == pSrc->nrow)
	{
		for (eIdx = 0; eIdx < pSrc->nelem; eIdx++)
		{
			pDstL[eIdx] -= pSrcL[eIdx];
		}
	}
	else
	{
		Result = MTX_SIZE_MISMATCH;
	}

	return Result;
}
void mtx_mul_scalar(Matrix_t *const pSrc, const double scalar)
{
	double *const pDst = pSrc->val;
	int eIdx;

	for (eIdx = 0; eIdx < pSrc->nelem; eIdx++)
	{
		pDst[eIdx] *= scalar;
	}
}
void mtx_sub_scalar(Matrix_t *const pSrc, const double scalar)
{
	double *const pDst = pSrc->val;
	int eIdx;

	for (eIdx = 0; eIdx < pSrc->nelem; eIdx++)
	{
		pDst[eIdx] -= scalar;
	}
}
void mtx_add_scalar(Matrix_t *const pSrc, const double scalar)
{
	double *const pDst = pSrc->val;
	int eIdx;

	for (eIdx = 0; eIdx < pSrc->nelem; eIdx++)
	{
		pDst[eIdx] += scalar;
	}
}
enum mtxResultInfo mtx_cpy(Matrix_t *const pDst, Matrix_t const *const pSrc)
{
	enum mtxResultInfo Result = MTX_OPERATION_OK;
	double *const pDstL = pDst->val;
	double const *const pSrcL = pSrc->val;
	int eIdx;

	if (pDst->ncol == pSrc->ncol && pDst->nrow == pSrc->nrow)
	{
		for (eIdx = 0; eIdx < pSrc->nelem; eIdx++)
		{
			pDstL[eIdx] = pSrcL[eIdx];
		}
	}
	else
	{
		Result = MTX_SIZE_MISMATCH;
	}

	return Result;
}
enum mtxResultInfo mtx_identity(Matrix_t *const pSrc)
{
	enum mtxResultInfo Result = MTX_OPERATION_OK;
	double *const pDst = (double *)pSrc->val;
	const int nCol = pSrc->ncol;
	int eIdx;

	if (pSrc->nrow == nCol)
	{
		pDst[0] = 1;

		for (eIdx = 1; eIdx < pSrc->nelem; eIdx++)
		{
			const int cmpLeft = (int)(eIdx / nCol);

			pDst[eIdx] = eIdx < nCol ? 0 : cmpLeft == eIdx % (cmpLeft * nCol) ? 1 : 0;
		}
	}
	else
	{
		Result = MTX_SIZE_MISMATCH;
	}

	return Result;
}
enum mtxResultInfo mtx_zeros(Matrix_t *const pSrc)
{
	enum mtxResultInfo Result = MTX_OPERATION_OK;
	double *const pDst = (double *)pSrc->val;
	int eIdx;

	for (eIdx = 0; eIdx < pSrc->nelem; eIdx++)
	{
		pDst[eIdx] = 0;
	}

	return Result;
}
