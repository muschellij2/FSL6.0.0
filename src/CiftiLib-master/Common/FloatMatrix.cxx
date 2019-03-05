/*LICENSE_START*/ 
/*
 *  Copyright (c) 2014, Washington University School of Medicine
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without modification,
 *  are permitted provided that the following conditions are met:
 *
 *  1. Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 *  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 *  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 *  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "CiftiAssert.h"
#include "CiftiException.h"
#include "FloatMatrix.h"
#include "MatrixFunctions.h"

using namespace cifti;
using namespace std;

bool FloatMatrix::checkDimensions() const
{
   uint64_t rows = m_matrix.size(), cols;
   if (rows == 0) return true;//treat it as fine for now
   cols = m_matrix[0].size();
   for (uint64_t i = 1; i < rows; ++i)
   {
      if (m_matrix[i].size() != cols)
      {
         return false;
      }
   }
   return true;
}

FloatMatrix::FloatMatrix(const vector<vector<float> >& matrixIn)
{
   m_matrix = matrixIn;
   CiftiAssert(checkDimensions());
}

FloatMatrix::FloatMatrix(const int64_t& rows, const int64_t& cols)
{
    resize(rows, cols, true);
}

bool FloatMatrix::operator!=(const FloatMatrix& right) const
{
   return !(*this == right);
}

FloatMatrix FloatMatrix::operator*(const FloatMatrix& right) const
{
   FloatMatrix ret;
   MatrixFunctions::multiply<float, float, float, double>(m_matrix, right.m_matrix, ret.m_matrix);
   return ret;
}

FloatMatrix& FloatMatrix::operator*=(const FloatMatrix& right)
{
   MatrixFunctions::multiply<float, float, float, double>(m_matrix, right.m_matrix, m_matrix);//would need a copy anyway, so let it make the copy internally
   return *this;
}

FloatMatrix FloatMatrix::concatHoriz(const FloatMatrix& right) const
{
   FloatMatrix ret;
   MatrixFunctions::horizCat(m_matrix, right.m_matrix, ret.m_matrix);
   return ret;
}

FloatMatrix FloatMatrix::concatVert(const FloatMatrix& bottom) const
{
   FloatMatrix ret;
   MatrixFunctions::vertCat(m_matrix, bottom.m_matrix, ret.m_matrix);
   return ret;
}

FloatMatrix FloatMatrix::getRange(const int64_t firstRow, const int64_t afterLastRow, const int64_t firstCol, const int64_t afterLastCol) const
{
   FloatMatrix ret;
   MatrixFunctions::getChunk(firstRow, afterLastRow, firstCol, afterLastCol, m_matrix, ret.m_matrix);
   return ret;
}

FloatMatrix FloatMatrix::identity(const int64_t rows)
{
   FloatMatrix ret;
   MatrixFunctions::identity(rows, ret.m_matrix);
   return ret;
}

FloatMatrix FloatMatrix::inverse() const
{
   FloatMatrix ret;
   MatrixFunctions::inverse(m_matrix, ret.m_matrix);
   return ret;
}

FloatMatrix& FloatMatrix::operator*=(const float& right)
{
   MatrixFunctions::multiply(m_matrix, right, m_matrix);//internally makes a copy
   return *this;
}

FloatMatrix FloatMatrix::operator+(const FloatMatrix& right) const
{
   FloatMatrix ret;
   MatrixFunctions::add(m_matrix, right.m_matrix, ret.m_matrix);
   return ret;
}

FloatMatrix& FloatMatrix::operator+=(const FloatMatrix& right)
{
   MatrixFunctions::add(m_matrix, right.m_matrix, m_matrix);
   return *this;
}

FloatMatrix& FloatMatrix::operator+=(const float& right)
{
   MatrixFunctions::add(m_matrix, right, m_matrix);
   return *this;
}

FloatMatrix FloatMatrix::operator-(const FloatMatrix& right) const
{
   FloatMatrix ret;
   MatrixFunctions::subtract(m_matrix, right.m_matrix, ret.m_matrix);
   return ret;
}

FloatMatrix& FloatMatrix::operator-=(const FloatMatrix& right)
{
   MatrixFunctions::subtract(m_matrix, right.m_matrix, m_matrix);
   return *this;
}

FloatMatrix& FloatMatrix::operator-=(const float& right)
{
   MatrixFunctions::add(m_matrix, (-right), m_matrix);
   return *this;
}

FloatMatrix& FloatMatrix::operator/=(const float& right)
{
   return ((*this) *= 1.0f / right);
}

bool FloatMatrix::operator==(const FloatMatrix& right) const
{
   if (this == &right)
   {
      return true;//short circuit true on pointer equivalence
   }
   int64_t i, j, rows = (int64_t)m_matrix.size(), cols;
   if (rows != (int64_t)right.m_matrix.size())
   {
      return false;
   }
   if (rows == 0)
   {
      return true;//don't try to get the second dimension
   }
   cols = (int64_t)m_matrix[0].size();
   if (cols != (int64_t)right.m_matrix[0].size())
   {
      return false;
   }
   for (i = 0; i < rows; ++i)
   {
      for (j = 0; j < cols; ++j)
      {
         if (m_matrix[i][j] != right.m_matrix[i][j])
         {
            return false;
         }
      }
   }
   return true;
}

void FloatMatrix::getDimensions(int64_t& rows, int64_t& cols) const
{
   rows = (int64_t)m_matrix.size();
   if (rows == 0)
   {
      cols = 0;
   } else {
      cols = (int64_t)m_matrix[0].size();
   }
}

FloatMatrixRowRef FloatMatrix::operator[](const int64_t& index)
{
   CiftiAssert(index > -1 && index < (int64_t)m_matrix.size());
   FloatMatrixRowRef ret(m_matrix[index]);
   return ret;
}

ConstFloatMatrixRowRef FloatMatrix::operator[](const int64_t& index) const
{
   CiftiAssert(index > -1 && index < (int64_t)m_matrix.size());
   ConstFloatMatrixRowRef ret(m_matrix[index]);
   return ret;
}

FloatMatrix FloatMatrix::reducedRowEchelon() const
{
   FloatMatrix ret(*this);
   MatrixFunctions::rref(ret.m_matrix);
   return ret;
}

void FloatMatrix::resize(const int64_t rows, const int64_t cols, const bool destructive)
{
   MatrixFunctions::resize(rows, cols, m_matrix, destructive);
}

FloatMatrix FloatMatrix::transpose() const
{
   FloatMatrix ret;
   MatrixFunctions::transpose(m_matrix, ret.m_matrix);
   return ret;
}

FloatMatrix FloatMatrix::zeros(const int64_t rows, const int64_t cols)
{
   FloatMatrix ret;
   MatrixFunctions::zeros(rows, cols, ret.m_matrix);
   return ret;
}

FloatMatrix FloatMatrix::ones(const int64_t rows, const int64_t cols)
{
   FloatMatrix ret;
   MatrixFunctions::ones(rows, cols, ret.m_matrix);
   return ret;
}

const vector<vector<float> >& FloatMatrix::getMatrix() const
{
   return m_matrix;
}

void FloatMatrix::getAffineVectors(Vector3D& xvec, Vector3D& yvec, Vector3D& zvec, Vector3D& offset) const
{
    if (m_matrix.size() < 3 || m_matrix.size() > 4 || m_matrix[0].size() != 4)
    {
        throw CiftiException("getAffineVectors called on incorrectly sized matrix");
    }
    xvec[0] = m_matrix[0][0]; xvec[1] = m_matrix[1][0]; xvec[2] = m_matrix[2][0];
    yvec[0] = m_matrix[0][1]; yvec[1] = m_matrix[1][1]; yvec[2] = m_matrix[2][1];
    zvec[0] = m_matrix[0][2]; zvec[1] = m_matrix[1][2]; zvec[2] = m_matrix[2][2];
    offset[0] = m_matrix[0][3]; offset[1] = m_matrix[1][3]; offset[2] = m_matrix[2][3];
}

FloatMatrix FloatMatrix::operator-() const
{
   int64_t rows, cols;
   getDimensions(rows, cols);
   FloatMatrix ret = zeros(rows, cols);
   ret -= *this;
   return ret;
}

FloatMatrixRowRef::FloatMatrixRowRef(vector<float>& therow) : m_row(therow)
{
}

FloatMatrixRowRef& FloatMatrixRowRef::operator=(const FloatMatrixRowRef& right)
{
   if (&m_row == &(right.m_row))
   {//just in case vector isn't smart enough to check self assignment
      return *this;
   }
   CiftiAssert(m_row.size() == right.m_row.size());//maybe this should be an exception, not an assertion?
   m_row = right.m_row;
   return *this;
}

FloatMatrixRowRef& FloatMatrixRowRef::operator=(const float& right)
{
   for (int64_t i = 0; i < (int64_t)m_row.size(); ++i)
   {
      m_row[i] = right;
   }
   return *this;
}

float& FloatMatrixRowRef::operator[](const int64_t& index)
{
   CiftiAssert(index > -1 && index < (int64_t)m_row.size());//instead of segfaulting, explicitly check in debug
   return m_row[index];
}

FloatMatrixRowRef::FloatMatrixRowRef(FloatMatrixRowRef& right) : m_row(right.m_row)
{
}

FloatMatrixRowRef& FloatMatrixRowRef::operator=(const ConstFloatMatrixRowRef& right)
{
   if (&m_row == &(right.m_row))
   {//just in case vector isn't smart enough to check self assignment
      return *this;
   }
   CiftiAssert(m_row.size() == right.m_row.size());
   m_row = right.m_row;
   return *this;
}

const float& ConstFloatMatrixRowRef::operator[](const int64_t& index)
{
   CiftiAssert(index > -1 && index < (int64_t)m_row.size());//instead of segfaulting, explicitly check in debug
   return m_row[index];
}

ConstFloatMatrixRowRef::ConstFloatMatrixRowRef(const ConstFloatMatrixRowRef& right) : m_row(right.m_row)
{
}

ConstFloatMatrixRowRef::ConstFloatMatrixRowRef(const vector<float>& therow) : m_row(therow)
{
}
