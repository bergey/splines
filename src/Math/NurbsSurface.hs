{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleContexts #-}

module Math.NurbsSurface where

import Data.List
import Data.Maybe
import Control.Applicative

import Data.VectorSpace
import Math.Spline.BSpline
import Math.Spline.Knots

data BSplineSurface s v = BSplineSurface  (Knots s) (Knots s) [[v]]
                        deriving Show

newtype NurbsSurface v = NurbsSurface (BSplineSurface (Scalar v) (Scalar v, v))

uKnots :: NurbsSurface v -> Knots (Scalar v)
uKnots (NurbsSurface (BSplineSurface k _ _)) = k

vKnots :: NurbsSurface v -> Knots (Scalar v)
vKnots (NurbsSurface (BSplineSurface _ k _)) = k

controlPoints :: NurbsSurface v -> [[(Scalar v, v)]]
controlPoints (NurbsSurface (BSplineSurface _ _ cps)) = cps

uDegree :: NurbsSurface v -> Int
uDegree (NurbsSurface (BSplineSurface k _ cps)) = m - n where
  m = numKnots k - 1
  n = length cps

vDegree :: NurbsSurface v -> Int
vDegree (NurbsSurface (BSplineSurface _ k cps)) = m - n where
  m = numKnots k - 1
  n = length $ head cps

unH :: (Fractional (Scalar v), VectorSpace v) => (Scalar v, v) -> v
unH (w,v) = recip w *^ v

toH :: VectorSpace v => Scalar v -> v -> (Scalar v, v)
toH w v = (w, w *^ v)

-- | @evalSurface n u v@ is the point on the surface n with
--   parametric coördinates u, v
evalSurface
  :: (Fractional (Scalar b), Ord (Scalar b), VectorSpace b,
      VectorSpace (Scalar b), Scalar (Scalar b) ~ Scalar b) =>
     NurbsSurface b -> Scalar b -> Scalar b -> Maybe b
evalSurface n u v = unH <$> evalSurface' n u v

evalSurface'
  :: (Fractional (Scalar v), Ord (Scalar v), VectorSpace (Scalar v),
      VectorSpace v, Scalar (Scalar v) ~ Scalar v) =>
     NurbsSurface v -> Scalar v -> Scalar v -> Maybe (Scalar v, v)
evalSurface' n u v | isNothing uspan = Nothing
                   | isNothing vspan = Nothing
                   | otherwise = Just . head . head $ rightMult uP (transpose [vFuns]) where
  uspan = findSpan (uKnots n) u
  vspan = findSpan (vKnots n) v
  uDeg  = uDegree n
  vDeg  = vDegree n
  uFuns = basisFuns' uDeg (uKnots n) u
  vFuns = basisFuns' vDeg (vKnots n) v
  leftMult = genMMult (*^) (^+^) zeroV
  rightMult = genMMult (^*) (^+^) zeroV
  rows = take (uDeg + 1) $ drop (fromJust uspan - uDeg) $ controlPoints n
  cps = map (take (vDeg + 1) . drop (fromJust vspan - vDeg)) rows
  uP    = leftMult [uFuns] cps

-- | surfaceGrid evaluates the NurbsSurface over a grid of points.
--   The grid is uniformly spaced (on each axis) in u, v, but not, in general,
--   in R3.
surfaceGrid ::
  (Enum (Scalar v), Fractional (Scalar v),
      Ord (Scalar v), VectorSpace v, VectorSpace (Scalar v),
      Scalar (Scalar v) ~ Scalar v) =>
     NurbsSurface v -- ^ surface to be evaluated
     -> Int           -- ^ number of points to evaluate on first (u) axis
     -> Int           -- ^ number of points to evaluate on second (v) axis
     -> [[v]]        -- ^ each inner list shares a value of u
surfaceGrid n uCt vCt = map f us where
  f u = mapMaybe (evalSurface n u) vs
  us = ctRange (uKnots n) (uDegree n) uCt
  vs = ctRange (vKnots n) (vDegree n) vCt
  ctRange ks p ct = case knotDomain ks p of
    Nothing       -> []
    Just (lo, hi) -> [lo, lo+(hi-lo)/(fromIntegral ct - 1)..hi]

-- | Generalized matrix matrix multiply
genMMult  :: (a -> b -> c)
        -> (c -> c -> c)
        -> c
        -> [[a]]
        -> [[b]]
        -> [[c]]
genMMult mul add c0 arr brr =
  [[foldr add c0 $ zipWith mul a b | b <- transpose brr] | a <- arr]
