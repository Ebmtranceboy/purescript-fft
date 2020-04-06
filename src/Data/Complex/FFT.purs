module Data.Complex.FFT where

import Prelude

import Control.Monad.Rec.Class (Step(..), tailRec)
import Control.Monad.ST (ST)
import Data.Array (reverse, replicate, mapWithIndex, length, sortBy)
import Data.Array.ST (STArray, push, pushAll, peek, poke, empty, run)
import Data.Complex (Cartesian(..), conjugate)
import Data.Foldable (foldr)
import Data.Int (toStringAs, binary, fromStringAs, toNumber, round)
import Data.Int.Bits (shl)
import Data.List (fromFoldable, toUnfoldable, transpose) as Array
import Data.Maybe (fromJust)
import Data.String (length) as String
import Data.String (toCodePointArray, singleton)
import Math (sqrt, cos, atan2, pi, log)
import PRNG ((!!))
import Partial.Unsafe (unsafePartial)

type Complex = Cartesian Number

type Index = Int
type PowerOfTwo = Int
type ExponentOfTwo = Int

initialSort :: Index -> ExponentOfTwo -> Index
initialSort i e2 =
  let s2' = toStringAs binary i
      s2 =  foldr (<>) "" (replicate (e2 - String.length s2') "0") <> s2'
      c2 = singleton <$> (reverse $ toCodePointArray s2)
  in unsafePartial $ fromJust $ fromStringAs binary $ foldr (<>) "" c2

data Direction = Forward | Backward

organize :: forall a. ExponentOfTwo -> Array Complex
              -> ST a (STArray a Complex)
organize b zs = tailRec go {idx: 0, arrM: empty}
  where go {idx, arrM} =
          if idx == shl 1 b
            then Done arrM
            else Loop { idx: idx+1
                      , arrM: do
                            arr <- arrM
                            _ <- push (zs !! (initialSort idx b)) arr
                            pure arr
                      }

trigTable :: forall a. PowerOfTwo -> ST a (STArray a Complex)
trigTable n = tailRec go { idx: 0
                         , arrM: do
                            a <- empty
                            _ <- push (Cartesian 1.0 0.0) a
                            pure a
                         }
  where go {idx, arrM} =
          let csn = cos (pi / toNumber n)
              sn = sqrt (1.0 - csn * csn)
            in if idx == n
                then Done arrM
                else Loop { idx: idx + 1
                          , arrM: do
                              arr <- arrM
                              Cartesian x y <- unsafePeek idx arr
                              _ <- push (Cartesian (x * csn - y * sn)
                                                   (y * csn + x * sn)) arr
                              pure arr
                          }

unsafePeek :: forall a. Int -> STArray a Complex -> ST a Complex
unsafePeek idx arr = do
  mz <- peek idx arr
  pure $ unsafePartial $ fromJust mz

process :: forall a. Direction -> ExponentOfTwo -> Array Complex
                     -> Array Complex -> ST a (STArray a Complex)
process dir b fs complexp = tailRec goi { i: 0
                                        , arrM: do
                                                   a <- empty
                                                   _ <- pushAll fs a
                                                   pure a
                                        }
  where gok j count step {k, arrM} =
          if k == step
            then Done arrM
            else Loop { k: k+1
                      , arrM: do
                                arr <- arrM
                                let idx = 2*j*step+k
                                f0 <- unsafePeek idx arr
                                fk <- unsafePeek (idx+step) arr
                                let cexp = complexp !! (count*k)
                                let f1 = fk * (case dir of
                                                  Forward -> conjugate cexp
                                                  Backward -> cexp)
                                _ <- poke idx (f0+f1) arr
                                _ <- poke (idx+step) (f0-f1) arr
                                pure arr
                      }
        goj count step {j, arrM} =
          if j == count
            then Done arrM
            else Loop {j: j+1, arrM: tailRec (gok j count step) {k: 0, arrM}}
        goi {i, arrM} =
          if i == b
            then Done arrM
            else
              let count = shl 1 (b-i-1)
                  step = shl 1 i
              in Loop {i: i+1, arrM: tailRec (goj count step) {j: 0, arrM}}

-- | U_i = fft Forward ((u_n))
-- | sr : sample rate
-- | Te : sampling period                               =>  Te = 1 / sr
-- | N : number of sample
-- | fe : frequency precision                           => fe = sr / (N - 1)
-- | T : signal duration
-- | u_n : sample at time t = n Te
-- | u_{N-1} : last sample at time t = T                => T = (N-1) Te = 1 / fe
-- | U_i : bin at frequency f = i fe                    => i  < (n-1) / 2
-- | U_{N/2-1} : last useful bin at frequency f = sr/2
fft :: Direction -> Array Complex -> Array Complex
fft dir zs =
  let n = toNumber $ length zs
      b = round $ log n / log 2.0
      n2 = (shl 1 b) `div` 2
      fs = run (organize b zs)
      complexp = run (trigTable n2)
   in (\z -> (_ / sqrt n) <$> z) <$> run (process dir b fs complexp)

transpose :: forall a. Array (Array a) -> Array (Array a)
transpose arr = Array.toUnfoldable $
  Array.toUnfoldable <$> (Array.transpose $
    Array.fromFoldable $ Array.fromFoldable <$> arr)

fft2 :: Direction -> Array (Array Complex) -> Array (Array Complex)
fft2 dir zss = transpose $ fft dir <$> (transpose $ fft dir <$> zss)

type Freq = Int
type Bin = { re :: Number
           , im :: Number
           , freq :: Freq
           , mag :: Number
           , phase :: Number
           }

bin :: Freq -> Complex -> Bin
bin k (Cartesian re im) =
  { re
  , im
  , freq: k
  , mag: sqrt $ re * re + im * im
  , phase: atan2 im re
  }

sortByMagnitude :: Array Complex -> Array Bin
sortByMagnitude zs = sortBy (\a b -> compare b.mag a.mag) $
    (\{index, value} -> bin index value) <$> (
        mapWithIndex (\index value -> {index, value}) zs)
