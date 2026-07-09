import sys
from collections import defaultdict, deque
import matplotlib.pyplot as plt

#DATA STRUCTURES

def make_polygon(vertices):
    area=0.0
    n=len(vertices)
    for i in range(n):
        x1,y1=vertices[i]
        x2,y2=vertices[(i + 1) % n]
        area +=x1*y2-x2*y1
    area /=2.0
    winding =+1 if area >1e-12 else (-1 if area<-1e-12 else 0)
    return {"vertices": vertices, "winding": winding}

def make_floor_plan(outer,holes):
    return {"outer": outer, "holes": holes}

def total_vertices(G):
    n=len(G["outer"]["vertices"])
    for h in G["holes"]:
        n+=len(h["vertices"])
    return n

#INPUT PARSING

def parse_input(raw):
    tokens=raw.split()
    pos=0
    def ri():
        nonlocal pos; v= int(tokens[pos]); pos +=1; return v
    def rf():
        nonlocal pos; v=float(tokens[pos]); pos+= 1; return v
    def ring(n): return [(rf(), rf()) for _ in range(n)]

    T = ri();cases=[]
    for _ in range(T):
        outer =make_polygon(ring(ri()))
        holes=[make_polygon(ring(ri())) for _ in range(ri())]
        cases.append(make_floor_plan(outer,holes))
    return cases



EPS=1e-9

def signed_area(verts):
    area=0.0; n=len(verts)
    for i in range(n):
        x1, y1 = verts[i]; x2, y2 = verts[(i+1) % n]
        area += x1*y2 - x2*y1
    return area / 2.0

def ensure_ccw(v): return v if signed_area(v) > 0 else list(reversed(v))
def ensure_cw(v):  return v if signed_area(v) < 0 else list(reversed(v))

def cross(o,a,b):
    return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])

def strict_in_triangle(p,a,b,c):
    #True ONLY if p is STRICTLY inside triangle (a,b,c)
    d1 = cross(a,b,p); d2=cross(b,c,p); d3=cross(c,a,p)
    return ((d1>EPS and d2>EPS and d3>EPS) or (d1<-EPS and d2<-EPS and d3<-EPS))

def point_in_triangle(p,a,b,c): 
    d1=cross(a,b,p); d2=cross(b,c,p); d3=cross(c,a,p)
    return not ((d1< 0 or d2 <0 or d3< 0) and (d1 > 0 or d2> 0 or d3 >0))

def segments_intersect_proper(p1,p2,p3,p4):
    #True if p1-p2 and p3-p4 cross strictly (not at endpoints)
    d1 = cross(p3,p4,p1); d2 = cross(p3,p4,p2)
    d3 =cross(p1,p2,p3); d4=cross(p1,p2,p4)
    return (((d1>EPS and d2<-EPS) or (d1 < -EPS and d2 > EPS)) and
            ((d3> EPS and d4 <-EPS) or (d3 < -EPS and d4 > EPS)))

def is_ear(prev,curr,nxt,ring):
    #True if curr is a valid ear (convex, empty interior)
    if cross(prev,curr,nxt) <=EPS: return False
    for v in ring:
        if v == prev or v == curr or v == nxt: continue
        if strict_in_triangle(v,prev,curr,nxt): return False
    return True


def is_visible_by_index(p,q_idx,polygon):
    q = polygon[q_idx]
    n = len(polygon)
    for i in range(n):
        a = polygon[i]
        b = polygon[(i+1)%n]
        # Skip edges sharing p as an endpoint
        if a == p or b == p:
            continue
        # Skip the two edges adjacent to the SPECIFIC occurrence q_idx
        if i == q_idx or i == (q_idx-1) % n:
            continue
        if segments_intersect_proper(p,q,a,b):
            return False
    return True




def ray_x_segment(M,a,b):
    """
    for "hole vertex on outer boundary":
      Use a half-open interval [minY, maxY) with a small epsilon so
      that a vertex exactly at M.y is counted for one edge only,
      preventing double-counting and the "no hit" case when M is
      collinear with an outer vertex.
    """
    ay, by = a[1],b[1]

    # Normalise so ay <= by
    if ay>by:
        a,b = b,a
        ay,by = by,ay

    # Half-open: include lower endpoint, exclude upper.
    # This guarantees exactly one edge is hit when M.y equals a vertex y.
    if not (ay-EPS <=M[1] <by+EPS):
        return False, None
    if abs(by-ay)<EPS:          # horizontal segment — skip
        return False, None

    t = (M[1]-ay)/(by-ay)
    x_int = a[0]+ t*(b[0]-a[0])

    # Strictly to the right, but allow a tiny epsilon so a touching
    # boundary is still detected (the "> M[0]" check was the bug
    # when M sat exactly on the outer edge).
    if x_int < M[0]-EPS:
        return False, None

    return True,x_int


def merge_holes(outer_verts,holes_verts):
    """
    Merge each hole into the outer ring via a bridge edge.

    edge cases:
      1. Half-open ray/segment interval → handles holes whose rightmost
         vertex lies exactly on an outer boundary edge .
      2. Reflex-vertex promotion → handles non-convex outers where the
         candidate endpoint is blocked by a reflex vertex of the outer ring.
      3. Index-based visibility check → avoids wrong-occurrence bugs when
         a bridge vertex appears twice in the merged ring.
      4. Holes processed rightmost-first so later bridges don't cross
         earlier ones.
    """
    merged = list(outer_verts)

    # Process holes in decreasing order of their rightmost x so that
    # earlier bridges don't obstruct later ones.
    sorted_holes = sorted(holes_verts,key=lambda h: max(p[0] for p in h),reverse=True)

    for hole in sorted_holes:
        n=len(merged)

        mi=max(range(len(hole)), key=lambda i: (hole[i][0], hole[i][1]))
        M =hole[mi]

        best_x    = float('inf')
        best_ei0  = None
        best_x_int = None

        for i in range(n):
            a=merged[i]
            b=merged[(i+1) % n]
            hit, x_int = ray_x_segment(M, a, b)
            if hit and x_int < best_x - EPS:
                best_x = x_int
                best_ei0 = i
                best_x_int = x_int

        if best_ei0 is None:
            # Fallback: scan again with a looser check (large epsilon)
            # This guards against all-horizontal / degenerate polygons.
            for i in range(n):
                a = merged[i]
                b = merged[(i+1) % n]
                ay, by = (a[1], b[1]) if a[1] <= b[1] else (b[1], a[1])
                a2, b2 = (a, b) if a[1] <= b[1] else (b, a)
                if ay-1e-6 <=M[1] <= by + 1e-6 and abs(by-ay) > 1e-9:
                    t = (M[1]-ay) / (by-ay)
                    x_int = a2[0] + t*(b2[0]-a2[0])
                    if x_int >= M[0] - 1e-6 and x_int < best_x - EPS:
                        best_x  = x_int
                        best_ei0  = i
                        best_x_int = x_int

        if best_ei0 is None:
            raise ValueError(
                f"No bridge edge found for hole vertex M={M}. "
                "Check that the hole is inside the outer boundary "
                "and that outer vertices are in CCW order.")

        ei0 = best_ei0
        ei1 = (best_ei0+1)%n
        p0, p1 = merged[ei0],merged[ei1]

        # Prefer the endpoint with the larger x (closer to M's ray).
        # If that endpoint is not visible from M (blocked by a reflex
        # vertex of the outer ring), try the other endpoint, then the
        # exact intersection point.
        #
        # ALSO do reflex-vertex promotion: if any outer vertex lies
        # inside the triangle (M, I, V) and is reflex, it becomes the
        # new V candidate — otherwise the bridge would cross it.

        I = (best_x_int,M[1])   # intersection point on the outer edge

        # Candidate V: endpoint with larger x
        if p0[0]>=p1[0]:
            primary_idx,secondary_idx = ei0,ei1
        else:
            primary_idx,secondary_idx = ei1,ei0

        # Reflex-vertex promotion inside triangle (M, I, V_cand)
        V_cand_idx = primary_idx
        V_cand     = merged[V_cand_idx]

        for j in range(n):
            q = merged[j]
            if j == V_cand_idx or j == ei0 or j == ei1:
                continue
            if q[0] <= M[0] - EPS:
                continue
            if not point_in_triangle(q, M, I, V_cand):
                continue
            #checking weather q is reflex vertex of the outer ring
            prev_j = merged[(j-1) % n]
            next_j = merged[(j+1) % n]
            if cross(prev_j, q, next_j) >= -EPS:
                continue    # convex — no threat
            # Promote: prefer the one which is closest to the +x direction
            dx_old = V_cand[0] - M[0]; dy_old = V_cand[1] - M[1]
            dx_new = q[0] - M[0];      dy_new = q[1] - M[1]
            len_old = (dx_old**2 + dy_old**2) ** 0.5 + EPS
            len_new = (dx_new**2 + dy_new**2) ** 0.5 + EPS
            cos_old = dx_old / len_old
            cos_new = dx_new / len_new
            if cos_new > cos_old + EPS:
                V_cand_idx = j
                V_cand     = q

        #Try the promoted candidate first, then fallback options
        V_idx = None
        for cand_idx in [V_cand_idx, secondary_idx]:
            if is_visible_by_index(M, cand_idx, merged):
                V_idx = cand_idx
                break

        if V_idx is None:
            #insert the intersection point as a new vertex
            insert_at = ei1 if ei1 != 0 else len(merged)
            merged.insert(insert_at, I)
            V_idx = insert_at

        V = merged[V_idx]
        #Result: ..., V, M, hole[mi], hole[mi+1], ..., hole[mi-1], M, V, ...
        hole_rotated = hole[mi:] + hole[:mi]    # starts at M

        merged = (merged[:V_idx + 1] +     
                  hole_rotated        +     
                  [M]                 +     
                  [V]                 +
                  merged[V_idx + 1:])

    return merged


#TRIANGULATION

def ear_clip(ring):
    triangles = []
    coords    = list(ring)
    n_init    = len(coords)

    if n_init < 3:
        return triangles, coords


    prev_v = [(i - 1) % n_init for i in range(n_init)]
    next_v = [(i + 1) % n_init for i in range(n_init)]
    active = set(range(n_init))

    def coord(i): return coords[i]

    def is_ear_v(ci):
        #O(|active|): containment check against every remaining vertex
        pi = prev_v[ci]; ni = next_v[ci]
        p, c, n = coord(pi), coord(ci), coord(ni)
        if cross(p, c, n) <= EPS:
            return False
        for idx in active:
            if idx == pi or idx == ci or idx == ni:
                continue
            if strict_in_triangle(coord(idx), p, c, n):
                return False
        return True

    def ear_area_v(ci):
        #Area proxy — no sqrt, used only for comparison.
        pi = prev_v[ci]; ni = next_v[ci]
        p, c, n = coord(pi), coord(ci), coord(ni)
        return abs(cross(p, c, n))
    
    ear_status = {i: is_ear_v(i) for i in range(n_init)}

    #Shared unlink helper
    def clip(ci):
        pi = prev_v[ci]; ni = next_v[ci]
        triangles.append((pi, ci, ni))
        # O(1) unlink
        active.discard(ci)
        next_v[pi] = ni
        prev_v[ni] = pi
        # O(n) re-evaluation — only 2 vertices, not all n
        ear_status[pi] = is_ear_v(pi)
        ear_status[ni] = is_ear_v(ni)

    while len(active) > 3:

        #smallest-area valid ear  O(n) scan 
        best_ear = None          # (area, vertex_index)
        for ci in active:
            if ear_status[ci]:
                area = ear_area_v(ci)
                if best_ear is None or area < best_ear[0]:
                    best_ear = (area, ci)

        if best_ear is not None:
            clip(best_ear[1])
            continue

        #collinear or degenerate vertex 
        ear_found = False
        for ci in list(active):        
            pi = prev_v[ci]; ni = next_v[ci]
            if abs(cross(coord(pi), coord(ci), coord(ni))) < 1e-9:
                clip(ci)
                ear_found = True
                break

        if ear_found:
            continue

        #force-clip most-convex vertex 
        best_ci = max(active,key=lambda ci: cross(coord(prev_v[ci]),coord(ci),coord(next_v[ci])))
        clip(best_ci)

    #Final triangle (exactly 3 vertices remain) 
    if len(active) == 3:
        a = next(iter(active))         
        b = next_v[a]; c = next_v[b]   
        triangles.append((a, b, c))

    return triangles, coords


def triangulate(G):
    outer_verts = ensure_ccw(list(G["outer"]["vertices"]))
    holes_verts = [ensure_cw(list(h["vertices"])) for h in G["holes"]]
    merged      = merge_holes(outer_verts, holes_verts)
    return ear_clip(merged)


#COLOURING

def three_colour(triangles, coords):
    if not triangles:
        return {}

    edge_to_tris = defaultdict(list)
    for tidx, (a, b, c) in enumerate(triangles):
        for u, v in [(a,b),(b,c),(c,a)]:
            edge_to_tris[(min(u,v), max(u,v))].append(tidx)

    colour_map = {}         
    visited    = set()
    all_idx    = set(range(len(triangles)))

    while all_idx - visited:
        start   = next(iter(all_idx - visited))
        a, b, c = triangles[start]

        free = list({0,1,2} - {colour_map.get(v) for v in [a,b,c]
                                if v in colour_map})
        for v in [a, b, c]:
            if v not in colour_map:
                colour_map[v] = free.pop(0)

        visited.add(start)
        queue = deque([start])

        while queue:
            tidx     = queue.popleft()
            tri      = triangles[tidx]
            used     = {colour_map[v] for v in tri if v in colour_map}

            for v in tri:
                if v not in colour_map:
                    remaining = list({0,1,2} - used)
                    colour_map[v] = remaining[0]
                    used.add(remaining[0])

            for u, v in [(tri[0],tri[1]),(tri[1],tri[2]),(tri[2],tri[0])]:
                key = (min(u,v), max(u,v))
                for nbr in edge_to_tris[key]:
                    if nbr not in visited:
                        visited.add(nbr)
                        queue.append(nbr)

    return colour_map 

def pick_guards(colour_map, coords):
    """
    Returns (guard_coord_list, colour_id).
    Translates the winning index-class back to coordinates for output/plot.
    """
    classes = {0: [], 1: [], 2: []}
    for idx, col in colour_map.items():
        classes[col].append(coords[idx])         
    best = min(classes, key=lambda c: len(classes[c]))
    return classes[best], best



def guard_bound(G): return (total_vertices(G)+2*len(G["holes"])) //3


#OUTPUT

def fmt_pt(p): return f"({p[0]:.4f},{p[1]:.4f})"

def print_results(case_num, triangles, k, guards):
    if case_num > 1: print()
    for tri in triangles:
        print(fmt_pt(tri[0]), fmt_pt(tri[1]), fmt_pt(tri[2]))
    print(k)
    for g in guards: print(fmt_pt(g))


#PLOTTING

COLOURS    = {0: "tomato",    1: "limegreen",  2: "royalblue"}
COL_LABELS = {0: "Class 0",   1: "Class 1",    2: "Class 2"}


def base(ax, G):
    outer = G["outer"]["vertices"]
    ox = [p[0] for p in outer] + [outer[0][0]]
    oy = [p[1] for p in outer] + [outer[0][1]]
    ax.fill(ox, oy, color="steelblue", alpha=0.10, zorder=1)
    ax.plot(ox, oy, "b-", lw=2, zorder=2, label="Outer boundary")
    for i, (px, py) in enumerate(outer):
        ax.scatter(px, py, color="blue", s=30, zorder=3)
        ax.annotate(f"p{i}", (px, py), xytext=(4, 4),
                    textcoords="offset points", fontsize=7, color="navy")
    for hi, hole in enumerate(G["holes"]):
        hpts = hole["vertices"]
        hx = [p[0] for p in hpts] + [hpts[0][0]]
        hy = [p[1] for p in hpts] + [hpts[0][1]]
        ax.fill(hx, hy, color="dimgray", alpha=0.70, zorder=4,
                label="Hole" if hi == 0 else "_")
        ax.plot(hx, hy, "k-", lw=1.5, zorder=5)
        for j, (px, py) in enumerate(hpts):
            ax.scatter(px, py, color="black", s=18, zorder=6)
            ax.annotate(f"h{hi}v{j}", (px, py), xytext=(4, 4),
                        textcoords="offset points", fontsize=6, color="black")


def tris(ax, triangles):
    for tri in triangles:
        tx = [tri[0][0], tri[1][0], tri[2][0], tri[0][0]]
        ty = [tri[0][1], tri[1][1], tri[2][1], tri[0][1]]
        ax.fill(tx, ty, color="cyan", alpha=0.07, zorder=2)
        ax.plot(tx, ty, color="deepskyblue", lw=0.9, ls="--", alpha=0.85, zorder=3)


def plot_gallery(G, triangles, colour_map, guards, guard_col_id, case_num):
    n = total_vertices(G); k = guard_bound(G)
    fig, axes = plt.subplots(1, 3, figsize=(21, 6))
    fig.suptitle(
        f"Test Case {case_num}   n={n}   k=⌊n+2h/3⌋={k}   "
        f"guards={len(guards)}  (colour class {guard_col_id})",
        fontsize=12, fontweight="bold")

    for ax, (title, show_tri) in zip(axes, [
            ("Input Polygon", False),
            ("Triangulation + 3-Colouring", True),
            ("Guard Placement", True)]):
        ax.set_aspect("equal"); ax.set_title(title, fontsize=11)
        ax.grid(True, ls="--", alpha=0.3)
        base(ax, G)
        if show_tri: tris(ax, triangles)

        if title == "Triangulation + 3-Colouring":
            seen = set()
            for vertex, col_id in colour_map.items():
                lbl = COL_LABELS[col_id] if col_id not in seen else "_"
                seen.add(col_id)
                ax.scatter(vertex[0], vertex[1], color=COLOURS[col_id],
                           s=90, zorder=8, label=lbl,
                           edgecolors="black", linewidths=0.5)
            ax.legend(loc="upper right", fontsize=7)

        elif title == "Guard Placement":
            for vertex, col_id in colour_map.items():
                ax.scatter(vertex[0], vertex[1], color=COLOURS[col_id],
                           s=50, zorder=6, alpha=0.30,
                           edgecolors="black", linewidths=0.3)
            if guards:
                gx = [g[0] for g in guards]; gy = [g[1] for g in guards]
                ax.scatter(gx, gy, c="red", s=220, marker="*", zorder=9,
                           label=f"Guards ({len(guards)})")
                for i, (xi, yi) in enumerate(zip(gx, gy)):
                    ax.annotate(f"G{i+1}", (xi, yi), xytext=(5, 5),
                                textcoords="offset points",
                                fontsize=8, color="red", fontweight="bold")
            ax.legend(loc="upper right", fontsize=8)
        else:
            ax.legend(loc="upper right", fontsize=8)

    plt.tight_layout()
    fname = f"gallery_case_{case_num}.png"
    plt.savefig(fname, dpi=150, bbox_inches="tight")
    print(f"[plot saved -> {fname}]", file=sys.stderr)
    plt.close(fig)



def main():
    raw = sys.stdin.read()
    test_cases = parse_input(raw)
    for case_num, G in enumerate(test_cases, 1):
        triangles, coords    = triangulate(G)      
        colour_map           = three_colour(triangles, coords)
        guards, guard_col    = pick_guards(colour_map, coords)
        k                    = len(guards)

        coord_tris = [(coords[a], coords[b], coords[c]) for a, b, c in triangles]

        coord_colour_map = {coords[i]: col for i, col in colour_map.items()}

        print_results(case_num, coord_tris, k, guards)
        plot_gallery(G, coord_tris, coord_colour_map, guards, guard_col, case_num)
if __name__ == "__main__":
    main()