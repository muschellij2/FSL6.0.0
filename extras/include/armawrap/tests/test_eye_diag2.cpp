#include "tests.hpp"

int main(int argc, char *argv[]) {

  IdentityMatrix i(10);
  DiagonalMatrix d(10);

  cout << "Identity matrix: " << endl;
  cout << "i.Minimum:              " << i.Minimum()                << endl;
  cout << "i.Maximum:              " << i.Maximum()                << endl;
  cout << "i.MinimumAbsoluteValue: " << i.MinimumAbsoluteValue()   << endl;
  cout << "i.MaximumAbsoluteValue: " << i.MaximumAbsoluteValue()   << endl;
  cout << "i.SumSquare:            " << i.SumSquare()              << endl;
  cout << "i.SumAbsoluteValue:     " << i.SumAbsoluteValue()       << endl;
  cout << "i.Sum:                  " << i.Sum()                    << endl;
  cout << "i.Trace:                " << i.Trace()                  << endl;
  cout << "i.Norm1:                " << i.Norm1()                  << endl;
  cout << "i.NormInfinity:         " << i.NormInfinity()           << endl;
  cout << "i.NormFrobenius:        " << i.NormFrobenius()          << endl;
  cout << "i.Determinant:          " << i.Determinant()            << endl;
  cout << "i.LogDeterminant:       " << i.LogDeterminant().Value() << endl;
  cout << "i.IsZero:               " << i.IsZero()                 << endl;


  cout << endl;
  cout << "Identity matrix = 5:" << endl;
  i = 5;
  cout << "i.Minimum:              " << i.Minimum()                << endl;
  cout << "i.Maximum:              " << i.Maximum()                << endl;
  cout << "i.MinimumAbsoluteValue: " << i.MinimumAbsoluteValue()   << endl;
  cout << "i.MaximumAbsoluteValue: " << i.MaximumAbsoluteValue()   << endl;
  cout << "i.SumSquare:            " << i.SumSquare()              << endl;
  cout << "i.SumAbsoluteValue:     " << i.SumAbsoluteValue()       << endl;
  cout << "i.Sum:                  " << i.Sum()                    << endl;
  cout << "i.Trace:                " << i.Trace()                  << endl;
  cout << "i.Norm1:                " << i.Norm1()                  << endl;
  cout << "i.NormInfinity:         " << i.NormInfinity()           << endl;
  cout << "i.NormFrobenius:        " << i.NormFrobenius()          << endl;
  cout << "i.Determinant:          " << i.Determinant()            << endl;
  cout << "i.LogDeterminant:       " << i.LogDeterminant().Value() << endl;
  cout << "i.IsZero:               " << i.IsZero()                 << endl;

  cout << endl; 
  cout << "Identity matrix = -5:" << endl;
  i = -5;
  cout << "i.Minimum:              " << i.Minimum()                << endl;
  cout << "i.Maximum:              " << i.Maximum()                << endl;
  cout << "i.MinimumAbsoluteValue: " << i.MinimumAbsoluteValue()   << endl;
  cout << "i.MaximumAbsoluteValue: " << i.MaximumAbsoluteValue()   << endl;
  cout << "i.SumSquare:            " << i.SumSquare()              << endl;
  cout << "i.SumAbsoluteValue:     " << i.SumAbsoluteValue()       << endl;
  cout << "i.Sum:                  " << i.Sum()                    << endl;
  cout << "i.Trace:                " << i.Trace()                  << endl;
  cout << "i.Norm1:                " << i.Norm1()                  << endl;
  cout << "i.NormInfinity:         " << i.NormInfinity()           << endl;
  cout << "i.NormFrobenius:        " << i.NormFrobenius()          << endl;
  cout << "i.Determinant:          " << i.Determinant()            << endl;
  cout << "i.LogDeterminant:       " << i.LogDeterminant().Value() << endl;
  cout << "i.IsZero:               " << i.IsZero()                 << endl;

  cout << endl; 
  cout << "Identity matrix = 0:" << endl;
  i = 0;
  cout << "i.Minimum:              " << i.Minimum()                << endl;
  cout << "i.Maximum:              " << i.Maximum()                << endl;
  cout << "i.MinimumAbsoluteValue: " << i.MinimumAbsoluteValue()   << endl;
  cout << "i.MaximumAbsoluteValue: " << i.MaximumAbsoluteValue()   << endl;
  cout << "i.SumSquare:            " << i.SumSquare()              << endl;
  cout << "i.SumAbsoluteValue:     " << i.SumAbsoluteValue()       << endl;
  cout << "i.Sum:                  " << i.Sum()                    << endl;
  cout << "i.Trace:                " << i.Trace()                  << endl;
  cout << "i.Norm1:                " << i.Norm1()                  << endl;
  cout << "i.NormInfinity:         " << i.NormInfinity()           << endl;
  cout << "i.NormFrobenius:        " << i.NormFrobenius()          << endl;
  cout << "i.Determinant:          " << i.Determinant()            << endl;
  cout << "i.LogDeterminant:       " << i.LogDeterminant().Value() << endl;
  cout << "i.IsZero:               " << i.IsZero()                 << endl;
  
  
  cout << endl;
  cout << "Diagonal matrix = [1, 10]:" << endl;
  d << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8 << 9 << 10;
  cout << "d.Minimum:              " << d.Minimum()                << endl;
  cout << "d.Maximum:              " << d.Maximum()                << endl;
  cout << "d.MinimumAbsoluteValue: " << d.MinimumAbsoluteValue()   << endl;
  cout << "d.MaximumAbsoluteValue: " << d.MaximumAbsoluteValue()   << endl;
  cout << "d.SumSquare:            " << d.SumSquare()              << endl;
  cout << "d.SumAbsoluteValue:     " << d.SumAbsoluteValue()       << endl;
  cout << "d.Sum:                  " << d.Sum()                    << endl;
  cout << "d.Trace:                " << d.Trace()                  << endl;
  cout << "d.Norm1:                " << d.Norm1()                  << endl;
  cout << "d.NormInfinity:         " << d.NormInfinity()           << endl;
  cout << "d.NormFrobenius:        " << d.NormFrobenius()          << endl;
  cout << "d.Determinant:          " << d.Determinant()            << endl;
  cout << "d.LogDeterminant:       " << d.LogDeterminant().Value() << endl;
  cout << "d.IsZero:               " << d.IsZero()                 << endl;

  cout << endl;
  cout << "Diagonal matrix = [-5, 4]:" << endl;
  d << -5 << -4 << -3 << -2 << -1 << 0 << 1 << 2 << 3 << 4;
  cout << "d.Minimum:              " << d.Minimum()                << endl;
  cout << "d.Maximum:              " << d.Maximum()                << endl;
  cout << "d.MinimumAbsoluteValue: " << d.MinimumAbsoluteValue()   << endl;
  cout << "d.MaximumAbsoluteValue: " << d.MaximumAbsoluteValue()   << endl;
  cout << "d.SumSquare:            " << d.SumSquare()              << endl;
  cout << "d.SumAbsoluteValue:     " << d.SumAbsoluteValue()       << endl;
  cout << "d.Sum:                  " << d.Sum()                    << endl;
  cout << "d.Trace:                " << d.Trace()                  << endl;
  cout << "d.Norm1:                " << d.Norm1()                  << endl;
  cout << "d.NormInfinity:         " << d.NormInfinity()           << endl;
  cout << "d.NormFrobenius:        " << d.NormFrobenius()          << endl;
  cout << "d.Determinant:          " << d.Determinant()            << endl;
  cout << "d.LogDeterminant:       " << d.LogDeterminant().Value() << endl;
  cout << "d.IsZero:               " << d.IsZero()                 << endl;


  cout << endl;
  cout << "Diagonal matrix = -10:" << endl;
  d = -10;
  cout << "d.Minimum:              " << d.Minimum()                << endl;
  cout << "d.Maximum:              " << d.Maximum()                << endl;
  cout << "d.MinimumAbsoluteValue: " << d.MinimumAbsoluteValue()   << endl;
  cout << "d.MaximumAbsoluteValue: " << d.MaximumAbsoluteValue()   << endl;
  cout << "d.SumSquare:            " << d.SumSquare()              << endl;
  cout << "d.SumAbsoluteValue:     " << d.SumAbsoluteValue()       << endl;
  cout << "d.Sum:                  " << d.Sum()                    << endl;
  cout << "d.Trace:                " << d.Trace()                  << endl;
  cout << "d.Norm1:                " << d.Norm1()                  << endl;
  cout << "d.NormInfinity:         " << d.NormInfinity()           << endl;
  cout << "d.NormFrobenius:        " << d.NormFrobenius()          << endl;
  cout << "d.Determinant:          " << d.Determinant()            << endl;
  cout << "d.LogDeterminant:       " << d.LogDeterminant().Value() << endl;
  cout << "d.IsZero:               " << d.IsZero()                 << endl;

  cout << endl;
  cout << "Diagonal matrix = 12:" << endl;
  d = 12;
  cout << "d.Minimum:              " << d.Minimum()                << endl;
  cout << "d.Maximum:              " << d.Maximum()                << endl;
  cout << "d.MinimumAbsoluteValue: " << d.MinimumAbsoluteValue()   << endl;
  cout << "d.MaximumAbsoluteValue: " << d.MaximumAbsoluteValue()   << endl;
  cout << "d.SumSquare:            " << d.SumSquare()              << endl;
  cout << "d.SumAbsoluteValue:     " << d.SumAbsoluteValue()       << endl;
  cout << "d.Sum:                  " << d.Sum()                    << endl;
  cout << "d.Trace:                " << d.Trace()                  << endl;
  cout << "d.Norm1:                " << d.Norm1()                  << endl;
  cout << "d.NormInfinity:         " << d.NormInfinity()           << endl;
  cout << "d.NormFrobenius:        " << d.NormFrobenius()          << endl;
  cout << "d.Determinant:          " << d.Determinant()            << endl;
  cout << "d.LogDeterminant:       " << d.LogDeterminant().Value() << endl;
  cout << "d.IsZero:               " << d.IsZero()                 << endl;
}
