# Intersection Notes

- [Why_intersections?](#why_intersections)
- [Lines_in_Two_Dimensions](#lines_in_two_dimensions)
- [More_Than_Two_Lines_in_Two_Dimensions](#more_than_two_lines_in_two_dimensions)
- [Planes_in_Three_Dimensions](#planes_in_three_dimensions)
- [Intersection_of_Planes_in_3D](#intersection_of_planes_in_3d)
- [Rules_for_Manipulating_Equations](#rules_for_manipulating_equations)
- [Solving_a_System_of_Linear_Equations](#solving_a_system_of_linear_equations)
- [Special_Cases_in_Gaussian_Elimination](#special_cases_in_gaussian_elimination)
-[Summary](#summary)

## Why_Intersections?

![1](intersection-pics/wi-1.png)
![2](intersection-pics/wi-2.png)
![3](intersection-pics/wi-3.png)

Equations come from observed or modeled relationships between real-world quantities.

![4](intersection-pics/wi-4.png)
![5](intersection-pics/wi-5.png)

Beta = 1 - Stock perfectly correlated to market
Beta = 0 - Stock completely uncorrelated to market
Beta = -1 - Stock perfectly negatively correlated to market

Beta(portfolio) is a weighted average of the stocks in the portfolio)
- therefore Wa + Wb = 1

To minimize market risk we want a Beta value for our portfolio to equal Zero

![6](intersection-pics/wi-6.png)

Finding the intersection between these 2 linear relationships shows us the weights that will give our portfolio a beta value of Zero and minimum market risk.

![7](intersection-pics/wi-7.png)

Solving systems of equations:
- Economics
- Biology
- Chemistry
- Physics
- Engineering
- Computer Science


## Lines_in_Two_Dimensions

![1](intersection-pics/l2d-1.png)
![2](intersection-pics/l2d-2.png)

A given line has infinite base points (any point on the line) and direction vectors (can multiply by non-zero scalar and will remain on same line)

![3](intersection-pics/l2d-3.png)

A more general equation to represent a line:

![4](intersection-pics/l2d-4.png)

Direction Vector

![5](intersection-pics/l2d-5.png)

Changing k shifts the line off the origin

![6](intersection-pics/l2d-6.png)

## More_Than_Two_Lines_in_Two_Dimensions

![1](intersection-pics/mt2l-1.png)

Common Intersection

![3](intersection-pics/mt2l-3.png)

More Common is individual intersections for each line

![2](intersection-pics/mt2l-2.png)

## Planes_in_Three_Dimensions

![1](intersection-pics/p3d-1.png)

Equation as the dot product of two vectors

![2](intersection-pics/p3d-2.png)

Changing k shifts the vector but does not change direction--just like in 2D.

![3](intersection-pics/p3d-3.png)

Vector [A, B, C] is the normal vector (orthogonal) to the plane.  If 2 planes have the same normal vector then they are parallel.

![4](intersection-pics/p3d-4.png)
![5](intersection-pics/p3d-5.png)
![6](intersection-pics/p3d-6.png)


## Intersection_of_Planes_in_3D

![1](intersection-pics/ip3d-1.png)

Two intersecting planes intersection will be a line.  The cross product of the two normal vectors will not be zero because they are not parallel and it will be orthogonal to both.

Pick an arbitrary point on intersection line (x0)

![2](intersection-pics/ip3d-2.png)
![3](intersection-pics/ip3d-3.png)

Consequences of this:

![4](intersection-pics/ip3d-4.png)

With introduction of 3rd plane

![5](intersection-pics/ip3d-51.png)
![6](intersection-pics/ip3d-6.png)
![7](intersection-pics/ip3d-7.png)
![8](intersection-pics/ip3d-8.png)


## Rules_for_Manipulating_Equations

With 3 planes the algebra gets crazy!  And will become less generalizable to higher dimensions.

![](intersection-pics/rfme-1.png)

Algorithm for solving system of linear system:

![](intersection-pics/rfme-2.png)
![](intersection-pics/rfme-3.png)
![](intersection-pics/rfme-4.png)

If we add the 2 equations the we can deduce that the solution is preserved and the intersection will remain the same.  It will keep all previous points and will not contain any new points.  

![](intersection-pics/rfme-5.png)

Can also add a multiple (or an arbitrary amount) of one equation to another.

![](intersection-pics/rfme-6.png)
![](intersection-pics/rfme-7.png)

## Solving_a_System_of_Linear_Equations

Build a procedure to reduce a system of equations to as simple a form as possible. Ideally, finding the unique solution to a system if it exists.

**Gaussian Elimination Algorithm**

Procedure:

![](intersection-pics/sse-1.png)

Use operations to reduce variables in system of equations so that there is one less variable in each equation.

Operations Reminder:
- Swap order of equations
- Multiply an equation by a non-zero number
- Add a multiple of an equation to another

![](intersection-pics/sse-2.png)

Notice the triangular form of the system.  The leading variable does not appear in subsequent equations.

![](intersection-pics/sse-3.png)

Note: This manipulated system does **not** represent **the same planes** as the original system **but** it does represent the **same solution set**.  We have changed the system such that the common intersections remain the same.

Once we have reached triangular form we clear upwards.  In this example we solve for z then y then x.

Solve for z:

![](intersection-pics/sse-4.png)

Solve for y:

![](intersection-pics/sse-5.png)

Solve for x to get common intersection point:

![](intersection-pics/sse-6.png)

Check to see that solution solves original system

![](intersection-pics/sse-7.png)

Run through process again with plot to demonstrate transformations occurring with each operation.

![](intersection-pics/sse-8.png)
![](intersection-pics/sse-9.png)
![](intersection-pics/sse-10.png)
![](intersection-pics/sse-11.png)
![](intersection-pics/sse-12.png)
![](intersection-pics/sse-13.png)
![](intersection-pics/sse-14.png)


**Practice**
![](intersection-pics/sse-15.png)


## Special_Cases_in_Gaussian_Elimination

**No Solution:**

![](intersection-pics/scge-1.png)
![](intersection-pics/scge-2.png)
![](intersection-pics/scge-3.png)
![](intersection-pics/scge-4.png)
![](intersection-pics/scge-5.png)

This contradiction tells us that the system has no solution.

**Multiple Solutions**

![](intersection-pics/scge-6.png)
![](intersection-pics/scge-7.png)
![](intersection-pics/scge-8.png)
![](intersection-pics/scge-10.png)

The system would have the same solution if we dropped the last equation all together.

To solve this we can **parametrize** the solution set (like we did for lines).

In the example above, once we remove the 3rd equation, see that we cannot clear terms downwards anymore.

Instead, we will identify a **pivot variable**--a variable that is a leading term when the system is in triangular form. In this case **x** and **y** are the leading variables (since they are in the front of the equations).

Since **z** is not a leading variable in any of the equations it is a **free variable**.  Therefore, it will become a parameter in the parametrization of a solution set.

![](intersection-pics/scge-11.png)
![](intersection-pics/scge-12.png)
![](intersection-pics/scge-13.png)
![](intersection-pics/scge-14.png)
![](intersection-pics/scge-15.png)

If we pick points for z we will find points for x and y that will solve the system from our parametrization.

**More Practice**

No Solution:
![](intersection-pics/gep-2.png)

Unique Solution:
![](intersection-pics/gep-3.png)

Infinite Solutions:
![](intersection-pics/gep-4.png)

Systems with No Solutions:
![](intersection-pics/gep-5.png)

Systems with One Solution:
![](intersection-pics/gep-6.png)

Systems with Many Solutions:
![](intersection-pics/gep-7.png)

**Summary of Characterization Results**
1. A system is inconsistent if we find 0=k for k non-zero during Gaussian Elimination.
2. Not enough to count the number of equations (usually) or look for 0=0 to determine if infinitely many solutions.
3. A consistent system has a unique solution if each variable is a pivot variable.
4. The number of free variables equals the dimension of the solution set.


## Summary

- Geometry and algebra are both important parts of solving systems of equations.
- Gaussian elimination exploits geometric insight to give us a procedure for finding an algebraic solution to a system.
- Next Step: Looking at systems of equations form a new perspective.  As statements about transformations of space.
