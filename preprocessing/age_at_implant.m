function days_old=age_at_implant(AnimalMetadata)
% age_at_implant: age in days at time of implant
%
% Ryan Harvey

    days_old=datenum(num2str(AnimalMetadata.Surgery.Date),'yyyymmdd')...
        -datenum(num2str(AnimalMetadata.Animal.DateOfBirth),'yyyymmdd');
end