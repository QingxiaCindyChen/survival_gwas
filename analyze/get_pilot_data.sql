-- select a.GRID, c.GENDER_EPIC, c.DOB, min(a.AGE_AT_EVENT) as AGE_AT_FIRST_EVENT, max(a.AGE_AT_EVENT) as AGE_AT_LAST_EVENT
-- from ICD_CODES a
-- inner join LB_EXOME b
-- on a.GRID = b.GRID
-- inner join SD_RECORD c
-- on a.GRID = c.GRID
-- where b.EURO = 1
-- group by a.GRID, c.GENDER_EPIC, c.DOB;

select a.GRID, c.GENDER_EPIC, c.DOB, min(a.ENTRY_DATE) as FIRST_ENTRY_DATE, max(a.ENTRY_DATE) as LAST_ENTRY_DATE
from ICD_CODES a
inner join LB_EXOME b
on a.GRID = b.GRID
inner join SD_RECORD c
on a.GRID = c.GRID
where b.EURO = 1
group by a.GRID, c.GENDER_EPIC, c.DOB;

# make one table per phenotype

# gout
select a.GRID, a.CODE, a.ENTRY_DATE
from ICD_CODES a
inner join LB_EXOME b
on a.GRID = b.GRID
where a.CODE in ('274', '274.0', '274.00', '274.01', '274.02', '274.03', '274.1', '274.10', '274.11', '274.19', '274.8', '274.81', '274.82', '274.89', '274.9')
and b.EURO = 1
order by a.GRID, a.ENTRY_DATE;

# ra
select a.GRID, a.CODE, a.ENTRY_DATE
from ICD_CODES a
inner join LB_EXOME b
on a.GRID = b.GRID
where a.CODE in ('714.0', '714.1', '714.2', '714.81')
and b.EURO = 1
order by a.GRID, a.ENTRY_DATE;

# ms
select a.GRID, a.CODE, a.ENTRY_DATE
from ICD_CODES a
inner join LB_EXOME b
on a.GRID = b.GRID
where a.CODE in ('340')
and b.EURO = 1
order by a.GRID, a.ENTRY_DATE;

# prostate cancer
select a.GRID, a.CODE, a.ENTRY_DATE
from ICD_CODES a
inner join LB_EXOME b
on a.GRID = b.GRID
where a.CODE in ('185', '233.4', 'V10.46')
and b.EURO = 1
order by a.GRID, a.ENTRY_DATE;

# a-fib
select a.GRID, a.CODE, a.ENTRY_DATE
from ICD_CODES a
inner join LB_EXOME b
on a.GRID = b.GRID
where a.CODE in ('427.31')
and b.EURO = 1
order by a.GRID, a.ENTRY_DATE;

# alzheimer's
select a.GRID, a.CODE, a.ENTRY_DATE
from ICD_CODES a
inner join LB_EXOME b
on a.GRID = b.GRID
where a.CODE in ('331.0', '331.00')
and b.EURO = 1
order by a.GRID, a.ENTRY_DATE;




