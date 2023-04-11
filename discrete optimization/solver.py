class ListSolver:

    def __init__(self, n_, C_):
        self.n = n_
        self.C = C_

    def solve(self, weights, price):
        """
        return max_p: Double, answer: List[Int]
        """
        n, C = self.n, self.C
        flag = 0
        prev = 1
        W = 2
        S = 3
        p = []
        w = []
        start = [-1, -1, 0, 0]
        lists = [start]
        all_lists = [lists]
        for i_ in range(n):
            #w_, p_ = map(int, input().split())
            # p.append(p_)
            # w.append(w_)
            tmp_lists = []
            new_lists = []
            for l in range(len(lists)):
                if lists[l][W] + weights[i_] <= C:
                    new_lists.append([1, l, lists[l][W] + weights[i_], lists[l][S] + price[i_]])
            i = 0
            j = 0
            # print(f'{i_+1})', lists, '|', new_lists)
            while i < len(lists) and j < len(new_lists):
                # print(lists[i], 'vs', new_lists[j])
                if (lists[i][W] <= new_lists[j][W]) and (lists[i][S] >= new_lists[j][S]):
                    # print('list')
                    #tmp_lists.append([0, i, lists[i][W], lists[i][S]])
                    j += 1
                    #i += 1
                elif (lists[i][W] >= new_lists[j][W]) and (lists[i][S] <= new_lists[j][S]):
                    # print('new_list')
                    #tmp_lists.append(new_lists[j])
                    #j += 1
                    i += 1
                else:
                    # print('noone')
                    if lists[i][W] <= new_lists[j][W]:
                        tmp_lists.append([0, i, lists[i][W], lists[i][S]])
                        i += 1
                    else:
                        tmp_lists.append(new_lists[j])
                        j += 1

            if i < len(lists):
                for k in range(i, len(lists)):
                    tmp_lists.append([0, k, lists[k][W], lists[k][S]])
                    # print(f'add i {k}')
            if j < len(new_lists):
                for k in range(j, len(new_lists)):
                    tmp_lists.append(new_lists[k])
                    # print(f'add j {k}')

            #print(f'{i_+1})', tmp_lists)
            lists = tmp_lists
            all_lists.append(lists)

        best = lists[-1]

        # for l in all_lists:
        #     print(l)
        result = []
        pos = -1
        for i in range(len(all_lists) - 1, 0, -1):
            if all_lists[i][pos][flag]:
                result.append(i-1)
            pos = all_lists[i][pos][prev]

        return best[S], result
