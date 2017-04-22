class PermWeightedSum:
    def __init__(self, length, sum_value):
        self.vect = [0] * length
        self.length = length
        self.sum_value = sum_value
        self.results = set()

    def get_weighted_sum(self):
        return sum([(i+1)*v for i, v in enumerate(self.vect)])

    def add_vect(self):
        ret = self.get_weighted_sum()
        if ret == len(self.vect) and sum(self.vect) == self.sum_value:
            self.results.add(tuple(self.vect))
        return ret

    def get(self, left_sum, idx):
        if left_sum < 0:
            return
        elif left_sum == 0:
            self.add_vect()
            return
        for k in range(len(self.vect)+1):
            self.vect[idx] = k
            if idx > 0:
                self.get(left_sum - (idx+1)*k, idx-1)
                current_sum = self.get_weighted_sum()
            else:
                current_sum = self.add_vect()
            self.vect[idx] = 0
            if current_sum > len(self.vect):
                break

    def __iter__(self):
        if len(self.results) == 0:
            self.get(len(self.vect), self.length-1)
        for el in self.results:
            yield el
