# ContactTracing
Simulating retrospective contact tracing of EMOD-HIV epidemics

###Decisions that we have been making###

When in the epidemic do we do the contact tracing?
     ~2010, or anywhere in the epidemic decline, basically
How long is the CT "campaign/activity"?
     Starting with 1 day, just to see how it works, but then expanding to 3 to 5 years
     Meaning within that 3 to 5 years everyone alive in that timeframe is eligible to be an index
Who is eligible to be an index individual?
     All HIV+ individuals
          All HIV+ individuals who are newly diagnosed
newly diagnosed, e.g. potentially before the campaign started? or
               actually diagnosed within the campaign timeframe
          All HIV+, newly diagnosed, in early infection? Or regardless of infection stage, as long as they are not virally suppressed (which they wouldn't be if they are
          newly diagnosed
What is the window for potential partners?
     e.g. "Name/find all sexual partners from the past 3 years"
     
     
###Additional considerations###

If someone is a contact of multiple incident HIV cases what should we do?

If we are tracing from every incident person we should just consider them the contact of the earliest diagnosed incident person (that is the first time they will show up).

If we are interested in not tracing contacts of every incident case or not all contacts get traced we need to keep them as a contact until they are traced. This requires me resolving the many to many mapping problem when trying to get the number of times a contact transmits before and/or after the tracing.