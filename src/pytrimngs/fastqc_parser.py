import re

class FastQC_Parser:
    @classmethod
    def parse_fastqc_data(cls, fastqc_string):
        modules = {}
        last_module = None
        mod = []
        for line in fastqc_string.split('\n'):
            line = line.rstrip()
            if 'NaN' in line: continue  
            fields = line.split("\t")
            if fields[0] == ">>END_MODULE":
                modules[last_module] = mod
                last_module = None
                mod = []
            elif fields[0] == "##FastQC":
                continue
            elif re.match(r'^>>', fields[0]):
                last_module = re.sub('>>', '', fields[0])
            else:
                mod.append(fields)
        processed_modules = {}
        processed_modules['general_stats'] = cls.parse_general_stats(modules['Basic Statistics'])
        processed_modules['quality_per_base'] = cls.parse_quality_per_base(modules['Per base sequence quality'])
        processed_modules['quality_per_sequence'] = cls.parse_quality_per_sequence(modules['Per sequence quality scores'])
        processed_modules['indeterminations_per_base'] = cls.parse_two_column_data(modules['Per base N content'])
        processed_modules['sequence_length_distribution'] = cls.parse_two_column_data(modules['Sequence Length Distribution'])
        #puts processed_modules.inspect
        return processed_modules

    @classmethod
    def parse_general_stats(cls, general_stats):
        stats = {}
        general_stats.pop(0) # Remove header
        for attrib, value in general_stats:
            if attrib == "Total Sequences":
                stats[attrib] = int(value)
            elif attrib == '%GC':
                stats[attrib] = float(value)
            elif attrib == "Sequence length":
                if '-' in value:
                    min_val, max_val = value.split('-')
                    stats['Read_max_length'] = int(max_val)
                    stats['Read_min_length'] = int(min_val)
                else:
                    stats['Read_max_length'] = int(value)        
                    stats['Read_min_length'] = int(value)
            else:
                stats[attrib] = value
        return stats

    @classmethod
    def parse_quality_per_base(cls, quality_data):
        quality_data.pop(0)
        for i, data in enumerate(quality_data):
            rec = [data.pop(0)]
            rec.extend([ float(d) for d in data ])
            quality_data[i] = rec
        return quality_data

    @classmethod
    def parse_quality_per_sequence(cls, quality_data):
        quality_data.pop(0)
        for i, data in enumerate(quality_data): quality_data[i] = [ float(d) for d in data ]
        return quality_data

    @classmethod
    def parse_two_column_data(cls, data):
        data.pop(0)
        new_data = [ [base, float(count)] for base, count in data ]
        return new_data

    @classmethod
    def get_mean(cls, data, col):
        total = 0
        count = 0
        for d in data:
            total += d[col]
            count += 1
        return total/count

    @classmethod
    def get_min(cls, data, col):
        nums = [ d[col] for d in data ]
        return min(nums)

    @classmethod
    def get_weighted_mean(cls, data, col_val, col_weigth):
        total_weigth = 0
        sum_product = 0
        for d in data:
            total_weigth += d[col_weigth]
            sum_product += d[col_weigth] * d[col_val]
        return sum_product/total_weigth

    @classmethod
    def get_weighted_mean_with_intervals(cls, data, col_val, col_weigth):
        total_weigth = 0
        sum_product = 0
        for d in data:
            total_weigth += d[col_weigth]
            summ = 0
            for i in [ int(v) for v in d[col_val].split('-') ]: summ += i
            sum_product += d[col_weigth] * summ
        return sum_product/total_weigth

    @classmethod
    def parse_distributions(cls, two_column_table):
        distribution = [ f"{row[0]}\t{row[1]}" for row in two_column_table ]
        distribution_string = ':'.join(distribution)
        return distribution_string
