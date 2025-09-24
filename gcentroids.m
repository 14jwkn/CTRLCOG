% Code adapted from: Statistics and Machine Learning Toolbox

function [centroids, counts] = gcentroids(X, index, clusts, dist)
%GCENTROIDS Centroids and counts stratified by group.
p = size(X,2);
num = length(clusts);
if isfloat(X)
    centroids = NaN(num,p,'like',X);
    counts = zeros(num,1,'like',X);
else %X can be logic
    centroids = NaN(num,p);
    counts = zeros(num,1);
end

for i = 1:num
    members = (index == clusts(i));
    if any(members)
       counts(i) = sum(members);
       switch dist
            case 'sqeuclidean'
                centroids(i,:) = sum(X(members,:),1) / counts(i);
            case 'cityblock'
                % Separate out sorted coords for points in i'th cluster,
                % and use to compute a fast median, component-wise
                Xsorted = sort(X(members,:),1);
                nn = floor(.5*counts(i));
                if mod(counts(i),2) == 0
                    centroids(i,:) = .5 * (Xsorted(nn,:) + Xsorted(nn+1,:));
                else
                    centroids(i,:) = Xsorted(nn+1,:);
                end
            case {'cosine','correlation'}
                centroids(i,:) = sum(X(members,:),1) / counts(i); % unnormalized
            case 'hamming'
                % Compute a fast median for binary data, component-wise
                centroids(i,:) = 2*sum(X(members,:), 1) > counts(i);
        end
    end
end
