from convexhull import convexHull
from scipy import stats
from scipy import special

def h_measure(true_class, scores, severity_ratio=None, threshold=0.5, level=[0.95]):

    n_scores = float(len(scores))
    sum_true = float(sum(true_class))
    n0 = float(n_scores - sum_true)
    pi0 = float(n0 / n_scores)
    pi1 = float(sum_true / n_scores)

    sr = severity_ratio
    if sr is None:
        sr = pi1 / pi0

    sc = list(scores)
    sc.sort()
    out_scores = get_score_distributions(true_class, scores, sum_true, n0)

    auc = auc_solver(out_scores['s0'], out_scores['f1'], out_scores['s1'])
    switched = False
    criterion = auc < 0.2
    if criterion:
        switched = True
        s = [1.0 - i for i in scores]
        out_scores = get_score_distributions(true_class, scores, sum_true, n0)

    f1 = out_scores['f1']
    f0 = out_scores['f0']
    s0 = out_scores['s0']
    s1 = out_scores['s1']
    s_len = out_scores['s_len']

    # Misclassification stats
    misclass_dict = misclass_counts(scores, true_class, threshold)
    auc = auc_solver(out_scores['s0'], out_scores['f1'], out_scores['s1'])
    gini = 2 * auc - 1.0
    ks = max([abs(i - j) for i, j in zip(f0, f1)])
    cost_parameter = sr / (1.0 + sr)
    mer = min([pi0 * (1.0 - j) + pi1 * k for j, k in zip(f0, f1)])
    mwl = 2 * min([cost_parameter * pi0 * (1.0 - j) + (1.0 - cost_parameter) * pi1 * k for j, k in zip(f0, f1)])

    sens_fixed = []
    spec_fixed = []
    look_up_auc(f0, [1.0 - j for j in f1])
    for i in range(0, len(level)):
        sens_fixed += [look_up_auc(f0, [1.0 - j for j in f1], x=level[i])]
    for i in range(0, len(level)):
        spec_fixed += [look_up_auc(f1, f0, x=1.0 - level[i])]
    chull_points = convex_hull(f0, f1)
    g0 = chull_points[0]
    g1 = chull_points[1]
    hc = len(g0)
    sg0 = [0.0] + diff(g0)
    sg1 = [0.0] + diff(g1)
    auch = sum([si*(gi-0.5*sgi) for si, gi, sgi in zip(sg0, g1, sg1)])
    s_class0 = []
    s_class1 = []
    for si, yi in zip(scores, true_class):
        if yi == 0:
            s_class0 += [si]
        elif yi == 1:
            s_class1 += [si]
    s_class0.sort()
    s_class1.sort()
    cost = [si for si in range(1, (hc+2), 1)]
    print cost
    b0 = [si for si in range(2, hc+2, 1)]
    b1 = [si for si in range(2, hc+2, 1)]
    # print cost
    # print b0
    # print b1
    if sr > 0:
        shape1 = 2
        shape2 = 1 + (shape1-1) * (1/sr)
    else:
        shape1 = pi1 + 1
        shape2 = pi0 + 1
    cost[0] = 0
    cost[hc] = 1

    b00 = stats.beta(shape1, shape2)
    test = special.btdtr(0, 1, 1)
###################################################
def convex_hull(a, b):
    fa = [1.0-ind for ind in a]
    fb = [1.0 - ind for ind in b]
    fc = [max(ind, jnd) for ind, jnd in zip(fa, fb)]
    points = [(ind, jnd) for ind, jnd in zip(fa, fc)]
    ch = convexHull(points)
    return [it[0] for it in ch], [it[1] for it in ch]



def look_up_auc(x_curve, y_curve, x=0):
    result = None
    if all([True if i >= 0.0 else False for i in diff(x_curve)]):
        ind = None
        while ind is None:
            for i in range(0, len(x_curve)):
                if x_curve[i] - x > 0:
                    ind = i
                    break
        x1 = float(x_curve[ind - 1])
        x2 = float(x_curve[ind])
        y1 = float(y_curve[ind - 1])
        y2 = float(y_curve[ind])
        if x2 - x1 > 1.0:
            pos = (x2 - x) / (x2 - x1)
            result = (1.0 - pos) * y1 + pos * y2
        else:
            result = y2
    return result


def diff(curve):
    result = [curve[1] - curve[0]]
    for i in range(2, len(curve)):
        result += [curve[i] - curve[i - 1]]
    return result


def auc_solver(s0, f1, s1):
    auc_sum = 0.0
    for a, b, c in zip(s0, f1, s1):
        auc_sum += a * (b - 0.5 * c)
    return 1.0 - auc_sum


def get_score_distributions(y, s, n1, n0):
    s1_combined = zip(s, y)
    s1_combined.sort()
    s1 = [0.0] + [i / n1 for j, i in s1_combined] + [0.0]
    y2 = [1.0 - i for i in y]
    s0_combined = zip(s, y2)
    s0_combined.sort()
    s0 = [0.0] + [i / n0 for j, i in s0_combined] + [0.0]

    s_len = len(s1)

    f1 = list(cum_sum(s1))
    f0 = list(cum_sum(s0))
    return {'f1': f1, 'f0': f0, 's1': s1, 's0': s0, 's_len': s_len}


def cum_sum(li):
    total = 0.0
    for item in li:
        total += item
        yield total


def misclass_counts(pred_class, true_class, threshold):
    tp = 0.0
    fp = 0.0
    tn = 0.0
    fn = 0.0
    scaled_predicted = [1.0 if i > threshold else 0.0 for i in pred_class]
    for p, t in zip(scaled_predicted, true_class):
        if p == 1 and t == 1:
            tp += 1.0
        elif p == 1 and t == 0:
            fp += 1.0
        elif p == 0 and t == 0:
            tn += 1.0
        elif p == 0 and t == 1:
            fn += 1.0

    er = (fp + fn) / (tp + fp + tn + fn)
    sens = tp / (tp + fn)
    spec = tn / (tn + fp)
    prec = tp / (tp + fp)
    recall = sens
    tpr = recall
    fpr = 1 - spec
    f = 2 / (1 / prec + 1 / sens)
    youden = sens + spec - 1
    return {'tp': tp, 'fp': fp, 'tn': tn, 'fn': fn, 'er': er,
            'sens': sens, 'spec': spec, 'prec': prec, 'recall': recall,
            'tpr': tpr, 'fpr': fpr, 'f': f, 'youden': youden}


true_class = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1]
score = [1.8, 0.9, 1.6, 0.5, 1, 0.1, 0.2, 2.6, -0.4, -0.1]
print h_measure(true_class, score)
