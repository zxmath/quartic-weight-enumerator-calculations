"""Smooth quartic counts in rational-point normal forms."""

from collections import defaultdict
from itertools import product
from sage.all import GF, PolynomialRing, ProjectiveSpace, QQ, matrix, vector


def geometrically_smooth(c, g, ell, d, e, f, h, i, j, k, ring):
    y = ring.gen()
    infinity = c*y**2 + g*y + ell
    if infinity.gcd(infinity.derivative()).gcd(y**3 + f*y + k).degree() > 0:
        return False
    a = c*y**2
    b = d + e*y + f*y**2 + g*y**3
    p = b - a**2/3
    r = h + i*y + j*y**2 + k*y**3 + ell*y**4 - a*b/3 + 2*a**3/27
    if not p:
        return r.gcd(r.derivative()).degree() == 0
    if p.gcd(r).gcd(r.derivative()).degree() > 0:
        return False
    discriminant = 4*p**3 + 27*r**2
    if not discriminant:
        return False
    common = discriminant.gcd(2*p*r.derivative() - 3*r*p.derivative())
    while common.degree() > 0:
        supported_on_p = common.gcd(p)
        if supported_on_p.degree() == 0:
            return False
        common = common // supported_on_p
    return True


def normal_forms(field, contact):
    # Coefficient weights after x -> b^contact*x, y -> b*y, and division by b^contact.
    weights = {
        2: (4, 3, 2, 6, 5, 4, 3, 2),
        3: (6, 4, 9, 7, 5, 3, 1),
        4: (8, 5, 2, 12, 9, 6),
    }[contact]
    scalars = [b for b in field if b]
    multipliers = [tuple(b**w for w in weights) for b in scalars]
    for coefficients in product(field, repeat=len(weights)):
        stabilizer = 0
        for multiplier in multipliers:
            image = tuple(a*b for a,b in zip(coefficients,multiplier))
            if image < coefficients:
                break
            stabilizer += image == coefficients
        else:
            orbit_size = len(scalars) // stabilizer
            yield normalized_coefficients(field, contact, coefficients), orbit_size


def normalized_coefficients(field, contact, coefficients):
    if contact == 2:
        d,e,f,h,i,j,k,ell = coefficients
        c,g = field(1),field(0)
    elif contact == 3:
        d,e,h,i,j,k,ell = coefficients
        c,g,f = field(0),field(1),field(0)
    else:
        d,e,f,h,i,j = coefficients
        c,g,ell,k = field(0),field(0),field(1),field(0)
    return c,g,ell,d,e,f,h,i,j,k


def evaluation_data(field):
    plane = ProjectiveSpace(field, 2, names=('x','y','z'))
    x,y,z = plane.gens()
    monomials = (x*z**3, y**2*z**2, y**3*z, y**4,
                 x**3*z, x**2*y*z, x*y**2*z, x**4, x**3*y, x**2*y**2, x*y**3)
    points = [tuple(point) for point in plane.rational_points()]
    values = matrix(field, [[m(*point) for m in monomials] for point in points])
    second = [(0,0),(0,1),(0,2),(1,1),(1,2),(2,2)]
    variables = (x,y,z)
    second_derivatives = [[m.derivative(variables[a]).derivative(variables[b])
                           for m in monomials] for a,b in second]
    derivatives = [matrix(field, [[d2m(*point) for d2m in row] for point in points])
                   for row in second_derivatives]
    return points,monomials,values,derivatives


def rational_point_profile(parameters, data, field):
    points,monomials,values,derivatives = data
    coefficients = vector(field, (1,)+parameters)
    on_curve = [j for j,value in enumerate(values*coefficients) if not value]
    hessian_values = [D*coefficients for D in derivatives]
    ordinary = 0
    for j in on_curve:
        xx,xy,xz,yy,yz,zz = [v[j] for v in hessian_values]
        if matrix(field, [[xx,xy,xz],[xy,yy,yz],[xz,yz,zz]]).det():
            ordinary += 1
    if ordinary:
        return len(on_curve),ordinary,0
    F = sum(a*m for a,m in zip(coefficients,monomials))
    gradients = [F.derivative(v) for v in F.parent().gens()]
    t = PolynomialRing(field, 't').gen()
    flexes = 0
    for j in on_curve:
        point = points[j]
        gx,gy,gz = [dF(*point) for dF in gradients]
        if point[0]:
            direction = (0,gz,-gy)
        elif point[1]:
            direction = (gz,0,-gx)
        else:
            direction = (gy,-gx,0)
        restriction = F(*(a+t*b for a,b in zip(point,direction)))
        if restriction[3]:
            flexes += 1
    return len(on_curve),0,flexes


def smooth_distribution(q, progress=False):
    field = GF(q)
    if field.characteristic() <= 3:
        raise ValueError('Characteristic must be greater than three')
    ring = PolynomialRing(field, 'y')
    data = evaluation_data(field)
    weights = defaultdict(QQ)
    candidates = {}
    smooth = {}
    for contact in (2,3,4):
        candidates[contact] = 0
        smooth[contact] = 0
        for parameters,orbit_size in normal_forms(field,contact):
            candidates[contact] += 1
            if not geometrically_smooth(*parameters,ring):
                continue
            smooth[contact] += 1
            n,ordinary,flexes = rational_point_profile(parameters,data,field)
            if contact == 2:
                marked = ordinary
            elif contact == 3 and ordinary == 0:
                marked = flexes
            elif contact == 4 and ordinary == flexes == 0:
                marked = n
            else:
                continue
            assert marked > 0
            weights[n] += QQ(orbit_size)/((q-1)*marked)
        if progress:
            print('Contact',contact,':',candidates[contact],'orbits,',smooth[contact],'smooth',flush=True)
    return dict(sorted(weights.items())),candidates,smooth
