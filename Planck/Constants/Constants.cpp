#include "Constants.hpp"
#include <cmath>

namespace Planck {
      const matrix pauliX = {
      {complex(0, 0), complex(1, 0)},
      {complex(1, 0), complex(0, 0)}
      };

      const matrix pauliY = {
      {complex(0, 0), complex(0, -1)},
      {complex(0, 1), complex(0, 0)}
      };

      const matrix pauliZ = {
      {complex(1, 0), complex(0, 0)},
      {complex(0, 0), complex(-1, 0)}
      };

      const matrix S = {
      {complex(1, 0), complex(0, 0)},
      {complex(0, 0), complex(0, 1)}
      };

      const matrix SInv = {
      {complex(1, 0), complex(0, 0)},
      {complex(0, 0), complex(0, -1)}
      };

      const matrix hadamard = {
      {complex(1 / std::sqrt(2), 0), complex(1 / std::sqrt(2), 0)},
      {complex(1 / std::sqrt(2), 0), complex(-1 / std::sqrt(2), 0)}
      };

      const matrix CNOT = {
      {complex(1, 0), complex(0, 0), complex(0, 0), complex(0, 0)},
      {complex(0, 0), complex(1, 0), complex(0, 0), complex(0, 0)},
      {complex(0, 0), complex(0, 0), complex(0, 0), complex(1, 0)},
      {complex(0, 0), complex(0, 0), complex(1, 0), complex(0, 0)}
      };

      const matrix identity = {
      {complex(1, 0), complex(0, 0)},
      {complex(0, 0), complex(1, 0)}
      };

      const matrix p0 = {
      {complex(1, 0), complex(0, 0)},
      {complex(0, 0), complex(0, 0)}
      };

      const matrix p1 = {
      {complex(0, 0), complex(0, 0)},
      {complex(0, 0), complex(1, 0)}
      };
}